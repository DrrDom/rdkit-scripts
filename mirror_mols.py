#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
from rdkit import Chem
from read_input import read_input


def main():
    parser = argparse.ArgumentParser(description='Mirror molecules in XY plane. Useful for generation of mirror '
                                                 'conformers or axial/planar enantiomers.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file with 3D structures.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output SDF file with structures mirrored in XY plane.')
    parser.add_argument('-s', '--suffix', metavar='STRING', required=False, default='_1',
                        help='Specify the suffix which will be appended to all molecule names. Specify empty string if '
                             'do not want to alter names. Default: _1')

    args = parser.parse_args()

    w = Chem.SDWriter(args.output)
    for mol, mol_name in read_input(args.input):
        conf = mol.GetConformer()
        for a in mol.GetAtoms():
            i = a.GetIdx()
            pos = list(conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, pos[:2] + [-pos[2]])
        if args.suffix:
            mol.SetProp('_Name', mol.GetProp('_Name') + args.suffix)
        w.write(mol)


if __name__ == '__main__':
    main()
