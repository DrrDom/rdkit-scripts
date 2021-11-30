#!/usr/bin/env python3

__author__ = 'NAME SURNAME'

import argparse
import sys
from read_input import read_input
from rdkit import Chem
from multiprocessing import Pool, cpu_count


def calc(input_fname, output_fname, ncpu, verbose):

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    # MAIN CODE BLOCK HERE
    # pool is optional and can be removed
    # unnecessary imports should be removed as well


def main():
    parser = argparse.ArgumentParser(description='INSERT DESCRIPTION HERE.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SMILES file.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    calc(args.input, args.output, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
