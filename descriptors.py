#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
from multiprocessing import Pool
from read_input import read_input
from rdkit.Chem.AtomPairs.Pairs import GetAtomPairFingerprint


def calc_ap(items):
    mol, mol_name = items
    fp = GetAtomPairFingerprint(mol).GetNonzeroElements()
    return mol_name, fp


def main():
    parser = argparse.ArgumentParser(description='Calculate RDKit descriptors.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='SMILES or SDFfile.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='text file with a computed descriptor matrix.')
    parser.add_argument('-d', '--descr', metavar='STRING', default='AP', required=False,
                        help='type of descriptors: AP (atom pairs). Default: AP.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of cores for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    funcs = {'AP': calc_ap}

    pool = Pool(args.ncpu)

    mol_names = []
    fps = []
    for i, (mol_name, fp) in enumerate(pool.imap(funcs[args.descr], read_input(args.input), chunksize=100), 1):
        mol_names.append(mol_name)
        fps.append(fp)
        if i % 1000 == 0:
            sys.stdout.write(f'\r{i} molecules were processed')

    d = pd.DataFrame(fps).fillna(0).astype(int)
    d.index = mol_names
    d.index.name = 'mol'
    d.to_csv(args.output, sep='\t')


if __name__ == '__main__':
    main()
