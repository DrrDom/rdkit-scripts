#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse


def filter_mols(input_fname, output_fname, names_fname):

    with open(names_fname)as f:
        names = {line.strip() for line in f}

    with open(input_fname) as f_in, open(output_fname, 'wt') as f_out:
        molstr = []
        for line in f_in:
            molstr.append(line)
            if line.strip() == '$$$$':
                if molstr[0].strip() in names:
                    f_out.writelines(molstr)
                molstr = []


def main():
    parser = argparse.ArgumentParser(description='Filter input molecules by names.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SDF file.')
    parser.add_argument('-n', '--names', metavar='FILENAME', required=True, type=str,
                        help='text file with molecule names to keep.')
    args = parser.parse_args()

    filter_mols(args.input, args.output, args.names)


if __name__ == '__main__':
    main()
