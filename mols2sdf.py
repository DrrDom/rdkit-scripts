#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 20-06-2018
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2018
# license         : 
#==============================================================================

__author__ = 'Pavel Polishchuk'

import os
import sys
import argparse
from rdkit import Chem


def main(input_fnames, output_fname):

    w = Chem.SDWriter(output_fname)

    for f in input_fnames:
        m = Chem.MolFromMolFile(f, sanitize=False, removeHs=False, strictParsing=False)
        if m:
            name = os.path.basename(f).rsplit('.', 1)[0]
            m.SetProp('Code', name)
            m.SetProp('_Name', name)
            w.write(m)
        else:
            sys.stderr.write('Molecule in %s cannot be read.' % f)

    w.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert multiple MOL files to SDF. No checks are performed.')
    parser.add_argument('-i', '--input', metavar='input.mol', required=True, nargs='*',
                        help='input MOL files.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output SDF file.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fnames = v
        if o == "output": output_fname = v

    main(input_fnames, output_fname)
