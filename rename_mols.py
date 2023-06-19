#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from collections import defaultdict, OrderedDict

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D

from read_input import read_input


def merge_confs(mols):

    if len(mols) == 0:
        return None
    if len(mols) == 1:
        return mols[0]

    main_mol = mols[0]
    for mol in mols[1:]:
        ids = main_mol.GetSubstructMatch(mol, useChirality=True)
        if not ids:
            sys.stderr.write(f'Molecule {mol.GetProp("_Name")} does not match main molecule {main_mol.GetProp("_Name")}. Skipped.\n')
            continue
        for c in mol.GetConformers():
            pos = c.GetPositions()
            for query_id, atom_id in enumerate(ids):
                x, y, z = pos[query_id,]
                c.SetAtomPosition(atom_id, Point3D(x, y, z))
            main_mol.AddConformer(c, assignId=True)

    return main_mol


def rename_molecules(input_fname, output_fname, names_fname, prefix):

    mols_dict = defaultdict(list)  # inchi: [mol, ...]
    mol_names = OrderedDict()  # inchi: new_name

    with open(names_fname, 'wt') as f_names:
        for mol, mol_name in read_input(input_fname):
            if mol.GetConformer().Is3D():
                Chem.AssignStereochemistryFrom3D(mol)
            inchi = Chem.inchi.MolToInchi(mol)
            if inchi not in mol_names:
                mol_names[inchi] = f'{prefix}{str(len(mol_names) + 1).zfill(8)}'
            mols_dict[inchi].append(Chem.RemoveHs(mol))
            f_names.write(f'{mol_names[inchi]}\t{mol_name}\n')

    w = Chem.SDWriter(output_fname)
    try:
        for inchi, mols in mols_dict.items():
            mol = merge_confs(mols)
            if mol:
                mol.SetProp('_Name', mol_names[inchi])
                for c in mol.GetConformers():
                    w.write(mol, confId=c.GetId())
    finally:
        w.close()


def main():
    parser = argparse.ArgumentParser(description='Identify identical structures (conformers) in an SDF file and '
                                                 'rename them identically in the same manner prefix + 8 sequential '
                                                 'digits. atoms in conformers will be reordered in a consistent manner '
                                                 'across all conformers of the same molecule.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SDF file. The order of molecules will change corresponding to new molecule names '
                             'to group identical molecules (conformers).')
    parser.add_argument('-n', '--names', metavar='FILENAME', required=True, type=str,
                        help='output text file with correspondence of new and old names. If some names occur multiple '
                             'times in input SDF, they will be treated separately and will result in multiple output '
                             'lines. The order of molecule names will be kept the same as in input SDF.')
    parser.add_argument('-p', '--prefix', metavar='STRING', required=False, type=str, default='MOL',
                        help='prefix to add to molecule names. Default: MOL.')
    args = parser.parse_args()

    rename_molecules(args.input, args.output, args.names, args.prefix)


if __name__ == '__main__':
    main()
