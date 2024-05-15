#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 08-06-2017
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2017
# license         : GPL3
#==============================================================================

import argparse
from rdkit import Chem


def main_params(input_sdf_fname, field_name, output_sdf_fname, insert_smi):

    w = Chem.SDWriter(output_sdf_fname)
    w.SetKekulize(False)
    i = 1

    for m in Chem.SDMolSupplier(input_sdf_fname, False, False, False):
        if m:
            if field_name and m.HasProp(field_name):
                m.SetProp('_Name', m.GetProp(field_name))
            elif insert_smi:
                tmp = Chem.Mol(m)
                Chem.SanitizeMol(tmp)
                m.SetProp('_Name', Chem.MolToSmiles(tmp, isomericSmiles=True))
            else:
                m.SetProp('_Name', 'MolID_%i' % i)
                i += 1
            w.write(m)

    w.close()


def main():

    parser = argparse.ArgumentParser(description='Insert mol titles to sdf file from a specified field or insert '
                                                 'canonical SMILES. If no options were selected compounds will be '
                                                 'enumerated sequentially with MolID_ prefix.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input sdf file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output sdf file.')
    parser.add_argument('-f', '--field_name', default=None,
                        help='field name in input sdf file which will be copied to molecule title.')
    parser.add_argument('--smi', action='store_true', default=False,
                        help='insert canonical SMILES in molecule title.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "output": output_fname = v
        if o == "field_name": field_name = v
        if o == "smi": insert_smi = v

    main_params(input_fname, field_name, output_fname, insert_smi)


if __name__ == '__main__':
    main()
