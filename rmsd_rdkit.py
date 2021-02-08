#!/usr/bin/env python2

__author__ = 'Pavel Polishchuk'


from rdkit import Chem
import argparse
import sys
import os
from read_input import read_pdbqt


def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)


def rmsd(mol, ref):
    match_indices = mol.GetSubstructMatches(ref, uniquify=False, useChirality=True)
    if not match_indices:
        return None
    else:
        ref_coord = get_coord(ref)
        min_rmsd = float('inf')
        for ids in match_indices:
            mol_coord = get_coord(mol, ids)
            s = 0
            for r, m in zip(ref_coord, mol_coord):
                s += ((r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2) ** 0.5
            s = s / len(ref_coord)
            if s < min_rmsd:
                min_rmsd = s
        return round(min_rmsd, 3)


def main_params(input_fnames, input_smi, output_fname, ref_name, refsmi):

    if ref_name.endswith('.mol2'):
        ref = Chem.MolFromMol2File(ref_name)
    elif ref_name.endswith('.pdbqt'):
        ref = read_pdbqt(ref_name, refsmi)[0]
    else:
        sys.stderr.write('Wrong format of the reference file. Only MOL2 and PDBQT files are allowed.')
        raise ValueError

    if output_fname is not None:
        sys.stdout = open(output_fname, 'wt')

    # read input smiles in dict: {name: smi}
    smis = dict()
    if input_smi is not None:
        with open(input_smi) as f:
            for line in f:
                values = line.strip().split()
                if len(values) >= 2:
                    smis[values[1]] = values[0]
                else:
                    sys.stderr.write(f'Line "{line}" in input smiles does not have two fields- SMILES and mol name. Skipped.')

    for in_fname in input_fnames:

        if in_fname.endswith('.mol2'):
            mols = [Chem.MolFromMol2File(in_fname)]
        elif in_fname.endswith('.pdbqt') or in_fname.endswith('.pdbqt_out'):
            mols = read_pdbqt(in_fname, smis[os.path.splitext(os.path.basename(in_fname))[0]])
        else:
            sys.stderr.write(f'Wrong format of the input file - {in_fname}. Only MOL2 and PDBQT files are allowed.')
            raise ValueError

        for i, mol in enumerate(mols, 1):
            if mol is None:
                print(f'{in_fname}\t{i}\tCannot read structure')
            else:
                mol_rmsd = rmsd(mol, ref)
                if mol_rmsd is not None:
                    print(f'{in_fname}\t{i}\t{mol_rmsd}')
                else:
                    print(f'{in_fname}\t{i}\tNo matches')


def main():

    parser = argparse.ArgumentParser(description='Calc RMSD between a reference molecule and docked poses.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='*',
                        help='input MOL2/PDBQT files to compare with a reference molecule.')
    parser.add_argument('--input_smi', metavar='FILENAME', required=False, default=None,
                        help='SMILES of input molecules if they are in PDBQT format. No header. Space- or '
                             'tab-separated. Molecule names should correspond to PDBQT file names.')
    parser.add_argument('-r', '--reference', metavar='FILENAME', required=True,
                        help='reference molecule (from X-ray complex structure) in MOL2/PDBQT format.')
    parser.add_argument('-s', '--refsmi', metavar='SMILES', required=False, default=None,
                        help='SMILES of the reference molecule. It requires only for PDBQT input to assign bond '
                             'orders.')
    parser.add_argument('-o', '--output', metavar='FILENAME',
                        help='output text file. If omitted output will be in stdout.')

    args = parser.parse_args()

    main_params(args.input, args.input_smi, args.output, args.reference, args.refsmi)


if __name__ == '__main__':
    main()
