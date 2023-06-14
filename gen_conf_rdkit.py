#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import os
import sys
import gzip
import argparse
import pickle
from itertools import combinations
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool, cpu_count
from read_input import read_input


def prep_input(fname, id_field_name, nconf, minimize, energy, rms, remove_h, seed):
    input_format = 'smi' if fname is None else None
    for mol, mol_name in read_input(fname, input_format=input_format, id_field_name=id_field_name):
        yield mol, mol_name, nconf, minimize, energy, rms, remove_h, seed


def map_gen_conf(args):
    return gen_confs(*args)


def remove_confs(mol, energy, rms, remove_h):

    if energy is None and rms is None:
        if remove_h:
            mol = Chem.RemoveHs(mol)
        return mol

    e = []
    for conf in mol.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        if ff is None:
            sys.stderr.write(f'Molecule {Chem.MolToSmiles(mol)} do not have all parameters MMFF\n')
            return
        e.append((conf.GetId(), ff.CalcEnergy()))
    e = sorted(e, key=lambda x: x[1])

    if not e:
        return

    remove_ids = []

    if energy is not None:
        keep_ids = [e[0][0]]
        for item in e[1:]:
            if item[1] - e[0][1] <= energy:
                keep_ids.append(item[0])
            else:
                remove_ids.append(item[0])
    else:
        keep_ids = [c.GetId() for c in mol.GetConformers()]

    if remove_h:
        mol = Chem.RemoveHs(mol)

    if rms is not None:
        rms_list = [(i1, i2, AllChem.GetConformerRMS(mol, i1, i2)) for i1, i2 in combinations(keep_ids, 2)]
        while any(item[2] < rms for item in rms_list):
            for item in rms_list:
                if item[2] < rms:
                    i = item[1]
                    remove_ids.append(i)
                    break
            rms_list = [item for item in rms_list if item[0] != i and item[1] != i]

    for cid in set(remove_ids):
        mol.RemoveConformer(cid)

    # reorder conformers by energy
    keep_ids = [i for i, en in e if i not in set(remove_ids)]  # ids is ordered because e is ordered
    for c in mol.GetConformers():
        c.SetId(c.GetId() + 100000)
    for c in mol.GetConformers():
        i = keep_ids.index(c.GetId() - 100000) + 1
        c.SetId(i)

    # save conf energy to mol
    en = [en for i, en in e if i not in set(remove_ids)]
    mol.SetProp('energy', ';'.join(f'{i} {v:.2f}' for i, v in enumerate(en, 1)))

    return mol


def gen_confs(mol, mol_name, nconf, minimize, energy, rms, remove_h, seed):
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", mol_name)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=nconf, maxAttempts=nconf*4, randomSeed=seed)
    if minimize:
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    mol = remove_confs(mol, energy, rms, remove_h)
    return mol_name, mol


def confs_to_string(mol, add_suffix):
    """
    Concatenate molblocks of all conformers is a single string. Add energy and suffix if available or requested.
    :param mol:
    :param add_suffix: bool
    :return:
    """
    string = ''
    e = {}
    if mol.HasProp('energy'):
        for item in mol.GetProp('energy').split(';'):
            conf_id, v = item.split(' ')
            e[int(conf_id)] = v
    for i, conf_id in enumerate(sorted(c.GetId() for c in mol.GetConformers()), 1):
        s = Chem.MolToMolBlock(mol, confId=conf_id)
        if add_suffix:
            header, ss = s.split('\n', 1)
            s = f'{header.strip()}_{i}\n' + ss
        if e:
            s += f'>  <energy>\n{e[conf_id]}\n\n'
        s += '$$$$\n'
        string += s
    return string


def main_params(in_fname, out_fname, id_field_name, nconf, minimize, energy, rms, remove_h, add_suffix,
                ncpu, seed, verbose):

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    output_file_type = None
    if out_fname is not None:

        if os.path.isfile(out_fname):
            os.remove(out_fname)

        if out_fname.lower().endswith('.sdf.gz'):
            writer = gzip.open(out_fname, 'a')
            output_file_type = 'sdf.gz'
        elif out_fname.lower().endswith('.sdf'):
            writer = open(out_fname, 'at')
            output_file_type = 'sdf'
        elif out_fname.lower().endswith('.pkl'):
            writer = open(out_fname, 'wb')
            output_file_type = 'pkl'
        else:
            raise Exception("Wrong output file format. Can be only SDF, SDF.GZ or PKL.")

    nprocess = min(cpu_count(), max(ncpu, 1))
    p = Pool(nprocess)

    try:
        for i, (mol_name, mol) in enumerate(p.imap_unordered(map_gen_conf,
                                                             prep_input(in_fname,
                                                                        id_field_name,
                                                                        nconf,
                                                                        minimize,
                                                                        energy,
                                                                        rms,
                                                                        remove_h,
                                                                        seed),
                                                             chunksize=1),
                                            1):
            if mol:
                if output_file_type == 'pkl':
                    pickle.dump((mol, mol_name), writer, -1)
                else:
                    mol.SetProp("_Name", mol_name)
                    string = confs_to_string(mol, add_suffix)
                    if string:   # wrong molecules (no valid conformers) will result in empty string
                        if out_fname is None:
                            sys.stdout.write(string)
                            sys.stdout.flush()
                        else:
                            writer.write(string.encode("ascii") if output_file_type == 'sdf.gz' else string)
            if verbose and i % 10 == 0:
                sys.stderr.write('\r%i molecules passed' % (i, ))
                sys.stderr.flush()

    finally:
        p.close()

    if out_fname is not None:
        writer.close()


def main():
    parser = argparse.ArgumentParser(description='Generate specified number of conformers using RDKit. '
                                                 'Conformers ids will be assigned according to conformer energy, '
                                                 'but the actual order of conformers on RDKit Mol object will not '
                                                 'be changed. Property _energy will be added to the Mol objects.')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=False, default=None,
                        help='input file with structures to generate conformers. Allowed formats SDF or SMILES. '
                             'if omitted STDIN will be used. STDIN takes only SMILES input (one or two columns).')
    parser.add_argument('-o', '--out', metavar='output.sdf', required=False, default=None,
                        help='output SDF file where conformers are stored. If extension will be SDF.GZ the output file '
                             'will be automatically gzipped. Alternatively for faster storage output can be stored '
                             'in a file with extension PKL. That is pickled storage of tuples (mol, mol_name). '
                             'If the output option will be omitted the output will be done to STDOUT in SDF format.')
    parser.add_argument('-d', '--id_field_name', metavar='field_name', default=None,
                        help='field name of compound ID in input SDF file. If omitted for sdf molecule titles '
                             'will be used or SMILES strings as names.')
    parser.add_argument('-n', '--nconf', metavar='conf_number', default=50,
                        help='number of generated conformers. Default: 50.')
    parser.add_argument('--no_minimize', action='store_true', default=False,
                        help='do not minimize energy of conformers using MMFF.')
    parser.add_argument('-e', '--energy_cutoff', metavar='VALUE', default=None,
                        help='conformers with energy difference from the lowest found one higher than the specified '
                             'value will be discarded. Default: None.')
    parser.add_argument('-r', '--rms', metavar='rms_threshold', default=None,
                        help='only conformers with RMS higher then threshold will be kept. '
                             'Default: None (keep all conformers).')
    parser.add_argument('-x', '--remove_h', action='store_true', default=False,
                        help='remove hydrogen atoms after generation of conformers and energy calculation (if it was '
                             'requested) but before calculation of rms and save to file.')
    parser.add_argument('-s', '--seed', metavar='random_seed', default=-1,
                        help='integer to init random number generator. Default: -1 (means no seed).')
    parser.add_argument('--suffix', action='store_true', default=False,
                        help='add a sequential conformer number to molecule name. Only applicable for '
                             'SDF and SDF.GZ outputs.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "id_field_name": id_field_name = v
        if o == "nconf": nconf = int(v)
        if o == "ncpu": ncpu = int(v)
        if o == "energy_cutoff": energy = float(v) if v is not None else None
        if o == "seed": seed = int(v)
        if o == "rms": rms = float(v) if v is not None else None
        if o == "verbose": verbose = v
        if o == "remove_h": remove_h = v
        if o == "no_minimize": minimize = not v
        if o == "suffix": add_suffix = v

    main_params(in_fname=in_fname,
                out_fname=out_fname,
                id_field_name=id_field_name,
                nconf=nconf,
                minimize=minimize,
                energy=energy,
                rms=rms,
                remove_h=remove_h,
                add_suffix=add_suffix,
                ncpu=ncpu,
                seed=seed,
                verbose=verbose)


if __name__ == '__main__':
    main()

