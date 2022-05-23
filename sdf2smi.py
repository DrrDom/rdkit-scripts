#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import argparse
from rdkit import Chem


def calc(input_fname, output_fname, field_name, extract_fields, sep):
    with open(output_fname, "wt") as f:
        for m in Chem.SDMolSupplier(input_fname):
            Chem.AssignStereochemistryFrom3D(m)
            if m:
                smi = Chem.MolToSmiles(m, isomericSmiles=True)
                if field_name:
                    n = m.GetProp(field_name)
                else:
                    n = m.GetProp("_Name")
                if extract_fields is not None:
                    fields = []
                    for f_name in extract_fields:
                        fields.append(m.GetProp(f_name))
                    f.write(smi + sep + n + sep + sep.join(fields) + '\n')
                else:
                    f.write(smi + sep + n + '\n')


def main():
    parser = argparse.ArgumentParser(description='Convert SDF to canonical isomeric SMILES with mol titles.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input SDF file.')
    parser.add_argument('-o', '--output', metavar='output.smi', required=True,
                        help='output SMILES file.')
    parser.add_argument('-f', '--field_name', metavar='FIELD_NAME', default=None,
                        help='sdf filed name which contains mol names. Optional. '
                             'If not specified mol titles will be used.')
    parser.add_argument('-e', '--extract_fields', metavar='FIELD_NAMES', default=None, nargs='*',
                        help='sdf fields to extract.')
    parser.add_argument('-s', '--sep', metavar='SEPARATOR', default='\t',
                        help='separator in output file. Default: tab.')


    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fname = v
        if o == "output": output_fname = v
        if o == "field_name": field_name = v
        if o == "extract_fields": extract_fields = v
        if o == "sep": sep = v

    calc(input_fname, output_fname, field_name, extract_fields, sep)


if __name__ == '__main__':
    main()
