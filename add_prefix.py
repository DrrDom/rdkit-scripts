#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import argparse


def main():
    parser = argparse.ArgumentParser(description='Add a prefix to molecule names in SDF file .')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output SDF file with updated mol names.')
    parser.add_argument('-p', '--prefix', metavar='STRING', default=True,
                        help='prefix to add to molecule names')

    args = parser.parse_args()

    with open(args.input) as f_in, open(args.output, 'wt') as f_out:
        molstr = []
        for line in f_in:
            molstr.append(line)
            if line.strip() == '$$$$':
                molstr[0] = args.prefix + molstr[0]
                f_out.writelines(molstr)
                molstr = []


if __name__ == '__main__':
    main()
