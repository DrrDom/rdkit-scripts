#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from read_input import read_input
from rdkit import Chem
from multiprocessing import Pool, cpu_count


def process_mol(mol, mol_name):
    try:
        frags = Chem.GetMolFrags(mol, asMols=True)
        output = ''.join(f'{Chem.MolToSmiles(frag)}\t{mol_name}_{i}\n' for i, frag in enumerate(frags))
        return output
    except:
        sys.stderr.write(f'molecule {mol_name} was omitted due to sanitization errors\n')
        return None


def process_mol_map(items):
    return process_mol(*items)


def main(input_fname, output_fname, ncpu, verbose):

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    with open(output_fname, 'wt') as f:
        for j, output_str in enumerate(pool.imap(process_mol_map, read_input(input_fname)), 1):
            if output_str:
                f.write(output_str)
            if verbose and j % 1000 == 0:
                sys.stderr.write(f'\r{j} molecules passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Save disconnected components of input molecules as '
                                                 'individual molecules with added suffix to the name.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SMILES file.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()
    main(args.input, args.output, args.ncpu, args.verbose)
