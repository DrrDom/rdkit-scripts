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
import sys
from rdkit import Chem


def main_params(in_fname, out_fname, field_name, force):

    w = Chem.SDWriter(out_fname)
    w.SetKekulize(False)

    for m in Chem.SDMolSupplier(in_fname, False, False, False):
        if m:
            if m.HasProp(field_name) and not force:
                sys.stderr.write(f'The field {field_name} exists. If you want to overwrite it please use --force '
                                 f'argument\n')
                exit()
            m.SetProp(field_name, m.GetProp('_Name'))
            w.write(m)

    w.close()


def main():
    parser = argparse.ArgumentParser(description='Extract field values from sdf-files.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output SDF file')
    parser.add_argument('-f', '--field_name', metavar='STRING', required=True,
                        help='field name where to insert molecule title. If the field exists, the script will be '
                             'interrupted.')
    parser.add_argument('--force', action='store_true', default=False,
                        help='force to overwrite existing field name.')


    args = parser.parse_args()

    main_params(args.input, args.output, args.field_name, args.force)


if __name__ == '__main__':
    main()
