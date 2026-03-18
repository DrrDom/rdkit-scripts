#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from read_input import read_input

# H, B, C, N, O, F, P, S, Cl, Br, I
organic_atoms = {1, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53}


def process_mol(items):
    mol, mol_name = items
    if mol is None:
        return None
    elm = set(a.GetAtomicNum() for a in mol.GetAtoms())
    is_organic = elm.issubset(organic_atoms)
    return f'{Chem.MolToSmiles(mol)}\t{mol_name}\n', is_organic


def calc(input_fname, organic_fname, inorganic_fname, ncpu, verbose):

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    w_organic = open(organic_fname, 'wt') if organic_fname else None
    w_inorganic = open(inorganic_fname, 'wt') if inorganic_fname else None

    n_organic = 0
    n_inorganic = 0

    for i, res in enumerate(pool.imap(process_mol, read_input(input_fname)), 1):
        if res is None:
            continue
        line, is_organic = res
        if is_organic:
            n_organic += 1
            if w_organic:
                w_organic.write(line)
        else:
            n_inorganic += 1
            if w_inorganic:
                w_inorganic.write(line)
        if verbose and i % 1000 == 0:
            sys.stderr.write(f'\r{i} structures processed')

    pool.terminate()

    if verbose:
        sys.stderr.write(f'\n{n_organic} organic, {n_inorganic} inorganic structures\n')

    if w_organic:
        w_organic.close()
    if w_inorganic:
        w_inorganic.close()


def main():
    parser = argparse.ArgumentParser(description='Split input structures into organic and inorganic sets. '
                                                 'Organic atoms are: H, B, C, N, O, F, P, S, Cl, Br, I.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output SMILES file for organic structures (only atoms from H, B, C, N, O, F, P, S, Cl, Br, I).')
    parser.add_argument('-d', '--inorganic', metavar='FILENAME', required=False, default=None,
                        help='output SMILES file for inorganic structures (contain atoms outside the organic set).')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    if args.output is None and args.inorganic is None:
        parser.error('at least one of --output or --inorganic must be specified.')

    calc(args.input, args.output, args.inorganic, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
