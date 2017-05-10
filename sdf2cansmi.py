#!/usr/bin/env python
# author          : Pavel
# date            : 06.03.15
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2015
# license         : GPL3
#==============================================================================

import argparse
from rdkit import Chem


def main_params(sdf_fname, smi_fname):
    f = Chem.SDMolSupplier(sdf_fname)
    with open(smi_fname, "wt") as smi:
        for mol in f:
            print(mol.GetProp("_Name"))
            smi.write(Chem.MolToSmiles(mol) + '\t' + mol.GetProp("_Name") + '\n')



def main():

    parser = argparse.ArgumentParser(description='Conversion of sdf to canonical SMILES with RDKit.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input sdf-file with compounds.')
    parser.add_argument('-o', '--output', metavar='output.smi', required=True,
                        help='output file with SMILES.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": sdf_fname = v
        if o == "output": smi_fname = v

    main_params(sdf_fname, smi_fname)


if __name__ == '__main__':
    main()
