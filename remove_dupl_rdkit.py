#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
from rdkit import Chem
from rdkit.Chem import inchi
from read_input import read_input


def get_inchi_key(mol, stereo):
    inchi_key = inchi.MolToInchiKey(mol)
    if not stereo:
        q = inchi_key.split('-')
        inchi_key = q[0] + '-' + q[2]  # remove middle part responsible for stereo and isotopes
    return inchi_key


def main_params(in_fname, dupl_fname, out_fname, ref_fname, stereo, ref_only):

    saver = Chem.SDWriter(out_fname)

    d = dict()  # key - inchi_key, value - list of names

    # read reference compounds
    if ref_fname is not None:
        for i, (mol, mol_name) in enumerate(read_input(ref_fname)):
            if mol is not None:
                inchi_key = get_inchi_key(mol, stereo)
                if inchi_key not in d.keys():
                    d[inchi_key] = [mol_name]

    for i, (mol, mol_name)in enumerate(read_input(in_fname)):
        if mol is not None:
            inchi_key = get_inchi_key(mol, stereo)
            if inchi_key not in d.keys():
                saver.write(mol)
                if not ref_only:
                    d[inchi_key] = [mol_name]
            else:
                d[inchi_key].append(mol_name)

    if dupl_fname is not None:
        with open(dupl_fname, 'wt') as f:
            for values in d.values():
                for v in values[1:]:
                    f.write(values[0] + '\t' + v + '\n')


def main():

    parser = argparse.ArgumentParser(description='Remove duplicates by inchi keys comparison.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input SDF, SMI, SDF.GZ or PKL file.')
    parser.add_argument('-r', '--reference', metavar='reference.sdf', required=False, default=None,
                        help='reference sdf file. If this file is specified than compounds from the input file '
                             'will be omitted if they are present in the reference file. If the input file contains '
                             'duplicated structures which are not in the reference file they will be also filtered. '
                             'To change this behavior set --reference_only argument. '
                             'Supported file formats are SDF, SMI, SDF.GZ or PKL.')
    parser.add_argument('-t', '--reference_only', action='store_true', default=False,
                        help='if set this flag only molecules duplicated with reference file will be omitted but not '
                             'duplicates with in input file.')
    parser.add_argument('-d', '--duplicates', metavar='duplicates.txt', required=False, default=None,
                        help='if specified names of removed duplicates will be stored in this file.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output SDF file with removed duplicates.')
    parser.add_argument('-s', '--stereo', action='store_true', default=False,
                        help='if set this flag stereoconfiguration will be considered during duplicate searching. '
                             'If set then isotopes will be considered as well.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fnames = v
        if o == "output": output_fname = v
        if o == "reference": ref_fname = v
        if o == "duplicates": dupl_fname = v
        if o == "stereo": stereo = v
        if o == "reference_only": ref_only = v

    main_params(input_fnames, dupl_fname, output_fname, ref_fname, stereo, ref_only)


if __name__ == '__main__':
    main()
