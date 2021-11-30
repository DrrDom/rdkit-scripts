#!/usr/bin/env python3

import argparse
from rdkit import Chem
from read_input import read_input


def calc(input_fname, output_fname):
    w = Chem.SDWriter(output_fname)
    for mol, mol_name in read_input(input_fname):
        mol.SetProp('_Name', mol_name)
        for conf in mol.GetConformers():
            w.write(mol, conf.GetId())
    w.close()


def main():
    parser = argparse.ArgumentParser(description='Convert multi-conformer PKL file to SDF.')
    parser.add_argument('-i', '--input', metavar='input.pkl', required=True,
                        help='input SDF file.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output pickle file.')

    args = parser.parse_args()
    calc(args.input, args.output)


if __name__ == '__main__':
    main()
