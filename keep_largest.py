#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
from rdkit import Chem
from read_input import read_input


def main_params(in_fname, out_fname):

    with open(out_fname, 'wt') as fo:

        for mol, mol_name in read_input(in_fname, sanitize=True):

            if mol is not None:

                frags = Chem.GetMolFrags(mol, asMols=True)
                max_hac = 0
                output_frag = None
                for frag in frags:
                    if frag.GetNumHeavyAtoms() > max_hac:
                        output_frag = frag
                fo.write(f'{Chem.MolToSmiles(output_frag, isomericSmiles=True)}\t{mol_name}\n')


def main():
    parser = argparse.ArgumentParser(description='Keep the largest fragment by the number of heavy atoms '
                                                 'in each compound record. If components have the same number a random '
                                                 'one will be selected.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SDF format.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output file in SMILES format.')

    args = parser.parse_args()

    if args.input == "/dev/stdin":
        args.input = None

    main_params(args.input, args.output)


if __name__ == '__main__':
    main()
