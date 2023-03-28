#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import argparse
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser(description='Convert SMILES to SDF with additional fields if they are named and '
                                                 'exist.')
    parser.add_argument('-i', '--input', metavar='input.smi', required=True,
                        help='input SMILES file.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output SDF file.')
    parser.add_argument('--header', action='store_true', default=False,
                        help='Set this argument if your input file contains additional fields which should be stored '
                             'to output SDF. By default, the script will try to guess whether the header exists and '
                             'if so the additional fields will be stored automatically.')
    parser.add_argument('-s', '--sep', metavar='SEPARATOR', default='\t',
                        help='separator in input file. Default: tab.')

    args = parser.parse_args()

    w = Chem.SDWriter(args.output)

    with open(args.input) as fin:
        line = fin.readline()
        if args.header or Chem.MolFromSmiles(line.strip().split(args.sep)[0]) is None:
            headers = line.strip().split(args.sep)[1:]
        else:
            headers = []
        if not headers:
            tmp = line.strip().split(args.sep)
            mol = Chem.MolFromSmiles(tmp[0])
            if mol:
                if len(tmp) > 1:
                    mol.SetProp('_Name', tmp[1])
                    mol.SetProp(headers[0], tmp[1])
                w.write(mol)
        for line in fin:
            tmp = line.strip().split(args.sep)
            mol = Chem.MolFromSmiles(tmp[0])
            if mol:
                if len(tmp) > 1:
                    mol.SetProp('_Name', tmp[1])
                    if headers:
                        for (field_name, value) in zip(headers, tmp[1:]):
                            mol.SetProp(field_name, value)
                w.write(mol)

    w.close()


if __name__ == '__main__':
    main()
