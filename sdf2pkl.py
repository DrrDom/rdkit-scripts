#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 06-07-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

import argparse
import pickle
from rdkit import Chem


def main_params(input_fname, output_fname):

    d = dict()
    for m in Chem.SDMolSupplier(input_fname):
        if m:
            mol_name = m.GetProp('_Name')
            if mol_name not in d:
                d[mol_name] = m
            else:
                d[mol_name].AddConformer(m.GetConformer(), assignId=True)

    with open(output_fname, 'wb') as f:
        for mol_name, mol in d.items():
            pickle.dump((mol, mol_name), f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert SDF to multi-conformer PKL file.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input SDF file.')
    parser.add_argument('-o', '--output', metavar='output.pkl', required=True,
                        help='output pickle file.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "output": output_fname = v

    main_params(input_fname, output_fname)
