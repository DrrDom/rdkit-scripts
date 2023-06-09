#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import random

from rdkit.Chem import rdMolDescriptors
from rdkit.SimDivFilters import rdSimDivPickers

from read_input import read_input


def pick(input_fname, output_fname, distance_threshold, seed):

    mols = list(read_input(input_fname))
    mols = sorted(mols, key=lambda x: x[1])
    random.seed(seed)
    random.shuffle(mols)
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 2048) for m, _ in mols]
    lp = rdSimDivPickers.LeaderPicker()
    picks = lp.LazyBitVectorPick(fps, len(fps), distance_threshold)
    with open(output_fname, 'wt') as f:
        f.write('\n'.join(sorted(mols[i][1] for i in picks)))


def main():
    parser = argparse.ArgumentParser(description='Selection of diverse subset of molecules using the sphere '
                                                 'exclusion algorithm.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='names of selected molecules.')
    parser.add_argument('-d', '--distance', metavar='NUMERIC', default=0.65, type=float,
                        help='distance threshold to separate selected molecules, Distance = 1 - similarity. '
                             'Default: 0.65 (recommended for morgan2 fingerprints).')
    parser.add_argument('-s', '--seed', metavar='INTEGER', default=0, type=int,
                        help='random seed to get reproducible output.')
    args = parser.parse_args()

    pick(args.input, args.output, args.distance, args.seed)


if __name__ == '__main__':
    main()
