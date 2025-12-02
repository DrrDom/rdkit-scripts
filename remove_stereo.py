#!/usr/bin/env python3

import argparse
import sys
from read_input import read_input
from rdkit import Chem
from multiprocessing import Pool, cpu_count


def remove_stereo(input_fname, output_fname, verbose):
    with open(output_fname, 'wt') as f:
        for i, (mol, mol_name) in enumerate(read_input(input_fname), 1):
            if mol:
                Chem.RemoveStereochemistry(mol)
                f.write(f'{Chem.MolToSmiles(mol, isomericSmiles=True)}\t{mol_name}\n')
            if verbose and i % 1000 == 0:
                sys.stderr.write(f'\rProcessed {i} molecules')
    if verbose:
        sys.stderr.write(f'\rProcessed {i} molecules\n')


def main():
    parser = argparse.ArgumentParser(description='Remove all stereo configurations from input molecules.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SMILES file.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    remove_stereo(args.input, args.output, args.verbose)


if __name__ == '__main__':
    main()
