#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 10-06-2020
# copyright       : Pavel Polishchuk 2020
# license         : GPL3
#==============================================================================

import re
import argparse

patt = re.compile('> {1,2}<(.*)>( +\([0-9]+\))?')


def strip_fields(molstr):
    out = []
    i = 0
    while i < len(molstr):
        out.append(molstr[i])
        if patt.fullmatch(molstr[i]):
            tmp = []
            while not patt.fullmatch(molstr[i + 1]) and molstr[i + 1] != '$$$$':
                tmp.append(molstr[i + 1])
                i += 1
            out.extend([line for i, line in enumerate(tmp) if line != '' or i == len(tmp) - 1])
        i += 1
    return out


def main_params(in_fname, out_fname):
    with open(in_fname) as ifs, open(out_fname, 'wt') as ofs:
        molstr = []
        for line in ifs:
            molstr.append(line.rstrip())
            if line.strip() == '$$$$':
                molstr = strip_fields(molstr)
                ofs.write('\n'.join(molstr) + '\n')
                molstr = []


def main():
    parser = argparse.ArgumentParser(description='Remove empty lines in multi-line field values.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output text file with values of specified fields')

    args = parser.parse_args()

    main_params(args.input, args.output)


if __name__ == '__main__':
    main()
