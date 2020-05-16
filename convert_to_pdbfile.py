#!/usr/bin/env python3
# author          : Alina Kutlushina
# date            : 15.05.2020
# license         : BSD-3
#==============================================================================

__author__ = 'alina'

import os
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from multiprocessing import Pool
from read_input import read_input


def read_molfiles(input, output_path):
    for mol, mol_name in read_input(input):
        yield mol, mol_name, output_path


def map_write_pdb(args):
    return write_pdb_file(*args)


def write_pdb_file(mol, mol_name, output_path):
    mol = Chem.AddHs(mol)
    if '2D' in Chem.MolToMolBlock(mol):
        AllChem.EmbedMolecule(mol)
    mol.SetProp('_Name', mol_name)
    Chem.MolToPDBFile(mol, os.path.join(output_path, mol_name + '.pdb'), flavor=4)


def convert_to_pdb(input, output_path, nprocess):
    """
    1. Add hydrogens
    2. Obtain a molecule coordinates if input molecules have 2D structures
    """

    p = Pool(nprocess)
    for _ in p.imap_unordered(map_write_pdb,
                              read_molfiles(input, output_path),
                              chunksize=10):
        continue


def entry_point():
    parser = argparse.ArgumentParser(description='Convert to a PDB file with the addition of hydrogens and '
                                                 'coordinates, if the initial structure has a 2d format',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input path to .smi or .sdf file')
    parser.add_argument('-o', '--output_path', type=str, required=True,
                        help='path to a folder where will store molecules in pdb format files')
    parser.add_argument('-c', '--ncpu', type=int, default=1,
                        help='number of cpu to use for calculation.')

    args = parser.parse_args()
    output_path = args.output_path
    os.makedirs(output_path, exist_ok=True)
    convert_to_pdb(input=os.path.abspath(args.input),
                   output_path=os.path.abspath(output_path),
                   nprocess=args.ncpu)


if __name__ == '__main__':
    entry_point()