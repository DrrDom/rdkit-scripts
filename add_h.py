#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import argparse
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser(description='Add hydrogens to input files in place.')
    parser.add_argument('-i', '--input', metavar='FILENAME', nargs='+', required=True,
                        help='input MOL files. Multiple files can be supplied.')

    args = parser.parse_args()

    for fname in args.input:

        if fname.endswith('.mol'):
            mol = Chem.MolFromMolFile(fname, removeHs=False)
            mol = Chem.AddHs(mol, addCoords=True)
            w = Chem.SDWriter(fname)
            w.write(mol)
            w.close()


if __name__ == '__main__':
    main()
