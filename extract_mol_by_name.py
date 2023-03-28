#!/usr/bin/env python

__author__ = 'Pavel Polishchuk'

import argparse
import re
import sys


def main():
    parser = argparse.ArgumentParser(description='Extract molecules by their names from SDF file.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output SDF file with selected molecules')
    parser.add_argument('-n', '--names', metavar='FILENAME', required=True,
                        help='input text file with molecule names to extract')
    parser.add_argument('-r', '--regex', metavar='REGEX_STRING', default=None,
                        help='regex to extract actual name from molecule title. Useful if partial matching of '
                             'molecular names are required. The first matching of the group will be used as a name. '
                             "Example: '^(.*)_[0-9]+$' will take ID1 from ID1_1 name to compare with supplied names.")

    args = parser.parse_args()

    names = set(line.strip() for line in open(args.names))

    with open(args.input) as fin:

        with open(args.output, 'wt') as fout:

            molstr = []

            for line in fin:
                molstr.append(line)
                if line.strip() == '$$$$':
                    if not args.regex:
                        name = molstr[0].strip()
                    else:
                        matches = re.findall(args.regex, molstr[0].strip())
                        if matches:
                            name = matches[0]
                        else:
                            sys.stderr.write(f'Molecule {molstr[0].strip()} does not match the given regex {args.regex}\n')
                            name = None
                    if name in names:
                        fout.writelines(molstr)
                    molstr = []


if __name__ == '__main__':
    main()
