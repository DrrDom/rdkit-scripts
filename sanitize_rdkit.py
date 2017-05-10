#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import sys
import argparse
from rdkit import Chem
from read_input import read_input


def main(in_fname, out_fname, error_fname, input_format, output_format):

    if out_fname is None:
        out_fname = '/dev/stdout'
    if error_fname is None:
        error_fname = '/dev/stderr'

    if output_format == 'sdf':
        wout = Chem.SDWriter(out_fname)
    else:
        print("Error output format.")
        exit()

    werr = Chem.SmilesWriter(error_fname, delimiter='\t', isomericSmiles=True)
    werr.SetProps(['Sanitization_error'])

    for mol, mol_name in read_input(in_fname, input_format=input_format, sanitize=False):

        if mol is not None:
            
            err = Chem.SanitizeMol(mol, catchErrors=True)
            
            if err:

                sys.stderr.write('Error %i sanitizing molecule %s\n' % (err, mol_name))
                sys.stderr.flush()
                mol.SetProp('Sanitization_error', str(err))
                werr.write(mol)

            else:

                wout.write(mol)

    werr.close()
    wout.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Removes Multi-component compounds and compounds with '
                                                 'non-organic atoms.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SDF format.')
    parser.add_argument('-f', '--input_format', metavar='sdf', required=False, default='sdf',
                        help='input file format. Default: sdf.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=False, default=None,
                        help='output file in SDF format. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-g', '--output_format', metavar='sdf', required=False, default='sdf',
                        help='output file format. Default: sdf.')
    parser.add_argument('-e', '--error_output', metavar='output.smi', required=False, default=None,
                        help='output SMILES file for compounds having sanitization errors. '
                             'If omitted output will be redirected to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "error_output": error_fname = v
        if o == "input_format": input_format = v
        if o == "output_format": output_format = v

    if in_fname == "/dev/stdin":
        in_fname = None

    main(in_fname, out_fname, error_fname, input_format, output_format)


