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


def prep_input(fname, id_field_name, nconf, energy, rms, seed):
    input_format = 'smi' if fname is None else None
    for mol, mol_name in read_input(fname, input_format=input_format, id_field_name=id_field_name):
        yield mol, mol_name, nconf, energy, rms, seed


def map_gen_conf(args):
    return gen_confs(*args)


def remove_confs(mol, energy, rms):

    if energy is None and rms is None:
        return

    e = []
    for conf in mol.GetConformers():
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf.GetId())
        if ff is None:
            print(Chem.MolToSmiles(mol))
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


def gen_confs(mol, mol_name, nconf, energy, rms, seed):
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", mol_name)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=nconf, maxAttempts=nconf*4, randomSeed=seed)
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    remove_confs(mol, energy, rms)
    return mol_name, mol


def main_params(in_fname, out_fname, id_field_name, nconf, energy, rms, ncpu, seed, verbose):

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
        for i, (mol_name, mol) in enumerate(p.imap_unordered(map_gen_conf, prep_input(in_fname, id_field_name, nconf, energy, rms, seed), chunksize=10), 1):
            if mol:
                if output_file_type == 'pkl':
                    pickle.dump((mol, mol_name), writer, -1)
                else:
                    mol.SetProp("_Name", mol_name)
                    e = dict()
                    if mol.HasProp('energy'):
                        for item in mol.GetProp('energy').split(';'):
                            conf_id, v = item.split(' ')
                            e[int(conf_id)] = v
                        string = "$$$$\n".join(Chem.MolToMolBlock(mol, confId=conf_id) + f'>  <energy>\n{e[conf_id]}\n\n' for conf_id in sorted(c.GetId() for c in mol.GetConformers()))
                    else:
                        string = "\n$$$$\n".join(Chem.MolToMolBlock(mol, confId=conf_id) for conf_id in sorted(c.GetId() for c in mol.GetConformers()))
                    if string:   # wrong molecules (no valid conformers) will result in empty string
                        string += "$$$$\n"
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
    parser.add_argument('-e', '--energy_cutoff', metavar='VALUE', default=None,
                        help='conformers with energy difference from the lowest found one higher than the specified '
                             'value will be discarded. Default: None.')
    parser.add_argument('-r', '--rms', metavar='rms_threshold', default=None,
                        help='only conformers with RMS higher then threshold will be kept. '
                             'Default: None (keep all conformers).')
    parser.add_argument('-s', '--seed', metavar='random_seed', default=-1,
                        help='integer to init random number generator. Default: -1 (means no seed).')
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

    main_params(in_fname=in_fname,
                out_fname=out_fname,
                id_field_name=id_field_name,
                nconf=nconf,
                energy=energy,
                rms=rms,
                ncpu=ncpu,
                seed=seed,
                verbose=verbose)


if __name__ == '__main__':
    main()

