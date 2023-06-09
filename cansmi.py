#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 26-07-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

import sys
import argparse
from read_input import read_input
from multiprocessing import Pool, cpu_count
from rdkit import Chem


def calc(items):
    # items is a tuple (mol, mol_name)
    if items[0].GetConformer().Is3D():
        Chem.AssignStereochemistryFrom3D(items[0])
    smi = Chem.MolToSmiles(items[0], isomericSmiles=True)
    return smi, items[1]


def main():

    parser = argparse.ArgumentParser(description='Conversion of input file to canonical SMILES with RDKit.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input file with compounds. Supported formats: SMILES (*.smi), '
                             'SDF (*.sdf, *.sdf.gz), Python pickled (*.pkl).')
    parser.add_argument('-o', '--output', metavar='output.smi', required=True,
                        help='output file with canonical SMILES.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of CPUs to use for computation.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    with open(args.output, 'wt') as f:
        if args.ncpu > 1:
            p = Pool(max(1, min(args.ncpu, cpu_count())))
            iterator = p.imap(calc, read_input(args.input))
        else:
            iterator = (calc(line) for line in read_input(args.input))
        for i, res in enumerate(iterator, 1):
            if res[0]:
                f.write('\t'.join(res) + '\n')
            if args.verbose and i % 10000 == 0:
                sys.stderr.write(f'\r{i} molecules passed')


if __name__ == '__main__':
    main()
