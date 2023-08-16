#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import numpy as np
import pandas as pd
import sys
from read_input import read_input
from rdkit import Chem


def process_file(input_fname, field_name, data_fname, output_fname, omit_missing):

    d = pd.read_csv(data_fname, sep="\t", header=0, index_col=0, dtype=str)
    # replace all values consisting of only whitespaces with NaN
    d = d.replace(r'^\s*$', np.nan, regex=True)

    w = Chem.SDWriter(output_fname)
    try:
        for mol in Chem.SDMolSupplier(input_fname, sanitize=False):
            mol_name = mol.GetProp(field_name) if field_name else mol.GetProp('_Name')
            prop_names = set(mol.GetPropNames())
            if omit_missing and mol_name not in d.index:
                continue
            if mol_name in d.index:
                for col_name in d.loc[mol_name].index:
                    if d.loc[mol_name][col_name] is not np.nan:
                        if col_name in prop_names:
                            sys.stderr.write(f'Property {col_name} for molecule {mol_name} was overwritten\n')
                        mol.SetProp(col_name, d.loc[mol_name][col_name])
            w.write(mol)
    finally:
        w.close()


def main():
    parser = argparse.ArgumentParser(description='Insert data from a text file to an SDF file as additional fields.'
                                                 ' Identification of matched records will be made by molecule name or '
                                                 'a given field containing a molecule name.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF file.')
    parser.add_argument('-f', '--field_name', metavar='STRING', required=False, type=str, default=None,
                        help='an optimal argument to specify the field name where molecule names are stored. Otherwise '
                             'molecule titles will be used.')
    parser.add_argument('-d', '--datafile', metavar='FILENAME', required=True, type=str,
                        help='input tab-separate text file with a header. The first column is identifier of a molecule '
                             '(molecule title)')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SDF file where molecules listed in the input text file get additional data '
                             'fields.')
    parser.add_argument('-x', '--omit_missing', action='store_true', default=False,
                        help='omit molecules missed in the input text file.')
    args = parser.parse_args()

    process_file(args.input, args.field_name, args.datafile, args.output, args.omit_missing)


if __name__ == '__main__':
    main()
