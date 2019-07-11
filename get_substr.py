#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 20-05-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

import sys
import argparse
from rdkit import Chem
from read_input import read_input
from multiprocessing import Pool, cpu_count


def match(args):
    # args - mol and mol_title
    mol, mol_title = args
    if mol.HasSubstructMatch(patt):
        return Chem.MolToSmiles(mol, isomericSmiles=True), mol_title
    else:
        return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Select compounds matching the specified substructure.')
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input SMILES/SDF file. SMILES file should contain no header.')
    parser.add_argument('-o', '--out', metavar='output.smi', required=True,
                        help='output file with SMILES and compound names.')
    parser.add_argument('-s', '--substr', metavar='SMARTS_STRING', required=True,
                        help='SMARTS string used for substructure matching.')
    parser.add_argument('-f', '--field_name', metavar='SDF_FIELD_NAME', default=None,
                        help='if sdf was passed as input the field name with molecule title can be optionally '
                             'specified.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1,
                        help='Number of CPU cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "substr": pattern = v
        if o == "field_name": field_name = v
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    patt = Chem.MolFromSmarts(pattern)

    if ncpu == 1:

        with open(out_fname, 'wt') as f:
            match_counter = 0
            for i, (mol, mol_title) in enumerate(read_input(in_fname, id_field_name=field_name)):
                if mol.HasSubstructMatch(patt):
                    f.write('%s\t%s\n' % (Chem.MolToSmiles(mol, isomericSmiles=True), mol_title))
                    f.flush()
                    match_counter += 1
                if verbose and (i + 1) % 100 == 0:
                    sys.stderr.write('\r%i molecules passed; matches: %i' % (i + 1, match_counter))
                    sys.stderr.flush()

    else:
        pool = Pool(max(1, min(ncpu, cpu_count())))
        with open(out_fname, 'wt') as f:
            match_counter = 0
            for i, res in enumerate(pool.imap(match, read_input(in_fname, id_field_name=field_name)), 1):
                if res:
                    smi, mol_title = res
                    f.write('%s\t%s\n' % (smi, mol_title))
                    f.flush()
                    match_counter += 1
                if verbose and i % 100 == 0:
                    sys.stderr.write('\r%i molecules passed; matches: %i' % (i, match_counter))
                    sys.stderr.flush()

