#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
from rdkit import Chem


def main_params(in_fname, dupl_fname, out_fname, ref_fname, stereo):

    saver = Chem.SDWriter(out_fname)

    d = dict()  # key - smiles, value - list of names

    # read reference compounds
    if ref_fname is not None:
        for i, mol in enumerate(Chem.SDMolSupplier(ref_fname)):
            if mol is not None:
                mol_name = mol.GetProp("_Name")
                if not mol_name:
                    mol_name = 'autoname_id_' + str(i)
                cansmi = Chem.MolToSmiles(mol, isomericSmiles=stereo)
                if cansmi not in d.keys():
                    d[cansmi] = [mol_name]

    for i, mol in enumerate(Chem.SDMolSupplier(in_fname)):
        if mol is not None:
            mol_name = mol.GetProp("_Name")
            if not mol_name:
                mol_name = 'autoname_id_' + str(i)
            cansmi = Chem.MolToSmiles(mol, isomericSmiles=stereo)
            if cansmi not in d.keys():
                saver.write(mol)
                d[cansmi] = [mol_name]
            else:
                d[cansmi].append(mol_name)

    if dupl_fname is not None:
        with open(dupl_fname, 'wt') as f:
            for values in d.values():
                for v in values[1:]:
                    f.write(values[0] + '\t' + v + '\n')


def main():

    parser = argparse.ArgumentParser(description='Remove duplicates by inchi comparison.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input sdf file.')
    parser.add_argument('-r', '--reference', metavar='reference.sdf', required=False, default=None,
                        help='reference sdf file. If this file is specified than compounds from the input file '
                             'will be omitted if they are present in the reference file. If the input file contains '
                             'duplicated structures which are not in the reference file they will be also filtered.')
    parser.add_argument('-d', '--duplicates', metavar='duplicates.txt', required=False, default=None,
                        help='if specified names of removed duplicates will be stored in this file.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output sdf file with removed duplicates.')
    parser.add_argument('-s', '--stereo', action='store_true', default=False,
                        help='if set this flag stereoconfiguration will be considered during duplicate searching.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fnames = v
        if o == "output": output_fname = v
        if o == "reference": ref_fname = v
        if o == "duplicates": dupl_fname = v
        if o == "stereo": stereo = v

    main_params(input_fnames, dupl_fname, output_fname, ref_fname, stereo)


if __name__ == '__main__':
    main()
