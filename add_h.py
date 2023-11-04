#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import argparse
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser(description='Add hydrogens to input files.')
    parser.add_argument('-i', '--input', metavar='FILENAME', nargs='+', required=True,
                        help='input MOL/SDF files. Multiple files can be supplied. MOL files are modified in place. '
                             'SDF files result in new SDF files with suffix _H in file name.')

    args = parser.parse_args()

    for fname in args.input:

        if fname.endswith('.mol'):
            mol = Chem.MolFromMolFile(fname, removeHs=False)
            mol = Chem.AddHs(mol, addCoords=True)
            Chem.MolToMolFile(mol, fname)

        if fname.endswith('.sdf'):
            mols = [Chem.AddHs(mol, addCoords=True) for mol in Chem.SDMolSupplier(fname, removeHs=False) if mol is not None]
            w = Chem.SDWriter(fname[:-4] + '_H.sdf')
            for mol in mols:
                w.write(mol)
            w.close()


if __name__ == '__main__':
    main()
