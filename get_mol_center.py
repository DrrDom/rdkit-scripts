#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import numpy as np
from read_input import read_input


def main():
    parser = argparse.ArgumentParser(description='Returns a geometrical center of input 3D molecules. '
                                                 'All explicit atoms will be considered.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file with 3D structures.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='text file with molecule name and three columns of x, y and z coordinates. No header.')

    args = parser.parse_args()

    with open(args.output, 'wt') as fout:
        for mol, mol_name in read_input(args.input):
            conf = mol.GetConformer()
            xyz = conf.GetPositions()
            res = np.mean(xyz, axis=0)
            fout.write(mol.GetProp('_Name') + '\t' + '\t'.join(map(str, res)) + '\n')


if __name__ == '__main__':
    main()
