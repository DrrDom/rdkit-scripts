#!/usr/bin/env python3

__author__ = 'NAME SURNAME'

import argparse
import sys
from read_input import read_input
from rdkit import Chem
from multiprocessing import Pool, cpu_count

pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")

def neutralize_atoms(mol):
    # https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def neutralize_item(item):
    mol, mol_name = item
    mol = neutralize_atoms(mol)
    output = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) + '\t' + mol_name + '\n'
    return output


def neutralize(input_fname, output_fname, ncpu, verbose):

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    with open(output_fname, "wt") as f:
        for i, line in enumerate(pool.imap(neutralize_item,
                              read_input(input_fname),
                              chunksize=1), 1):
            if line:
                f.write(line)
            if verbose and i % 1000 == 0:
                sys.stderr.write(f'\rProcessed {i}')
    if verbose:
        sys.stderr.write(f'\n')


def main():
    parser = argparse.ArgumentParser(description='Neutralize input structures.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SMILES file.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    neutralize(args.input, args.output, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()



#
#
#
# #!/usr/bin/env python
# # coding: utf-8
#
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit import RDLogger
# RDLogger.DisableLog('rdApp.*')
# from rdkit.Chem.MolStandardize import rdMolStandardize
#
# import os.path as osp
# import numpy as np
# import pandas as pd
# import sys
# import os
# import datetime
# from multiprocessing import Pool, cpu_count
# from functools import partial
#
#
# root = osp.abspath(osp.dirname(__file__))
#
#
# def neutralize_atoms(mol):
#     # https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
#     pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
#     at_matches = mol.GetSubstructMatches(pattern)
#     at_matches_list = [y[0] for y in at_matches]
#     if len(at_matches_list) > 0:
#         for at_idx in at_matches_list:
#             atom = mol.GetAtomWithIdx(at_idx)
#             chg = atom.GetFormalCharge()
#             hcount = atom.GetTotalNumHs()
#             atom.SetFormalCharge(0)
#             atom.SetNumExplicitHs(hcount - chg)
#             atom.UpdatePropertyCache()
#     return mol
#
#
#
# def predict_neutralized_mol(item):
#     print(f'{datetime.datetime.now()}; PID: {os.getpid()}; {item[1]} start')
#     try:
#         m, mol_name = item
#
#         mol = deepcopy(m)
#         mol_n = neutralize_atoms(Chem.RemoveHs(mol))
#         smi = Chem.MolToSmiles(mol_n, isomericSmiles=True)
#         smi_n = Chem.CanonSmiles(smi)
#         sys.stdout.flush()
#
#     except Exception as e:
#         print(e)
#         return None
#     return smi_n + '\t' + mol_name + '\n'
#
# def read_input(input_fname, sanitize=True):
#     with open(input_fname) as f:
#         for line in f:
#             tmp = line.strip().split()
#             mol = Chem.MolFromSmiles(tmp[0], sanitize=sanitize)
#             if mol is not None:
#                 if len(tmp) > 1:
#                     mol_title = tmp[1]
#                 else:
#                     mol_title = __get_smi_as_molname(mol)
#                 mol.SetProp('_Name', mol_title)
#                 print(f'read {mol_title}')
#                 sys.stdout.flush()
#                 yield mol, mol_title
#             else:
#                 print(f'read failed')
#                 sys.stdout.flush()
#
#
# def main(input_fname, output_fname):
#
#     pool = Pool(cpu_count())
#
#     with open(output_fname, "w") as f:
#         for line in pool.imap(partial(predict_neutralized_mol),
#                               read_input(input_fname),
#                               chunksize=1):
#             if line:
#                 f.write(line)
#                 f.flush()
#
#
# if __name__=="__main__":
#     input_fname = sys.argv[1]
#     output_fname = sys.argv[2]
#     main(input_fname, output_fname)
#
