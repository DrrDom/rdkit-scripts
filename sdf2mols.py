#!/usr/bin/env python3

import argparse
import os


def main():

    parser = argparse.ArgumentParser(description='Split SDF to MOL files with molecule names.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file with molecule names.')
    parser.add_argument('-o', '--output', metavar='DIRNAME', required=True,
                        help='dir name where to store output MOL files.')

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    with open(args.input) as f_in:
        molstr = []
        for line in f_in:
            if line.strip() == '$$$$':
                mol_name = molstr[0].strip()
                with open(os.path.join(args.output, mol_name + '.mol'), 'wt') as f_out:
                    for l in molstr:
                        f_out.write(l)
                        if l.strip() == 'M  END':
                            break
                molstr = []
            else:
                molstr.append(line)


if __name__ == '__main__':
    main()
