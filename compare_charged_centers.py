#!/usr/bin/env python3

import argparse
import sys
import traceback
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from itertools import combinations

from neutralize import neutralize_atoms


def get_atoms_within_radius(mol, center, radius):
    visited = {center}
    frontier = {center}
    for _ in range(radius):
        new = set()
        for a in frontier:
            for nbr in mol.GetAtomWithIdx(a).GetNeighbors():
                new.add(nbr.GetIdx())
        visited |= new
        frontier = new
    return list(visited)


def atom_centered_smiles(mol, atom_idx, radius):
    """
    Return SMILES of an atom and atoms away with up to a given number of bonds. The center atoms is the first in SMILES string
    :param mol:
    :param atom_idx:
    :param radius:
    :return:
    """
    mol.GetAtomWithIdx(atom_idx).SetProp('center', 'true')
    atoms = get_atoms_within_radius(mol, atom_idx, radius)
    bonds = []
    for i, j in combinations(atoms, 2):
        b = mol.GetBondBetweenAtoms(i, j)
        if b:
            bonds.append(b.GetIdx())
    submol = Chem.PathToSubmol(mol, bonds)
    center_idx = None
    for a in submol.GetAtoms():
        if a.HasProp('center') and a.GetProp('center') == 'true':
            center_idx = a.GetIdx()
            break
    return Chem.MolToSmiles(submol, rootedAtAtom=center_idx)


def compare_smiles(smi1, smi2):

    output = []

    try:
        mol1 = Chem.MolFromSmiles(smi1)
        mol2 = Chem.MolFromSmiles(smi2)

        data1 = []
        for atom in mol1.GetAtoms():
            data1.append(atom.GetFormalCharge())

        data2 = []
        for atom in mol2.GetAtoms():
            data2.append(atom.GetFormalCharge())

        mol1 = neutralize_atoms(mol1)
        mol2 = neutralize_atoms(mol2)

        matches = mol1.GetSubstructMatches(mol2)

        min_s = 1000
        min_match = None
        for match in matches:
            s = 0
            for j, i in enumerate(match):
                s = s + (data1[i] - data2[j])
            if abs(s) < min_s:
                min_s = s
                min_match = match

        radius = 2

        if min_match:
            data1 = [data1[i] for i in min_match]
            for i, j in enumerate(min_match):
                if data1[j] != data2[i]:
                    if data1[j] == 0:
                        output.append((atom_centered_smiles(mol1, j, radius), data1[j], data2[i]))
                    else:
                        output.append((atom_centered_smiles(mol2, i, radius), data1[j], data2[i]))
                elif data1[j] != 0:
                    output.append((atom_centered_smiles(mol1, j, radius), data1[j], data2[i]))
    except Exception as e:
        traceback.format_exc()

    return output


def process_line(line):
    # line should have SMILES, id, protonated SMILES 1, protonated SMILES 2
    output = []
    line = line.strip().split()
    # print(line)
    items = compare_smiles(line[2], line[3])
    if items:
        for item in items:
            output.append(line[0] + '\t' + line[1] + '\t' + '\t'.join(map(str, item)) + '\n')
    return output


def main():
    parser = argparse.ArgumentParser(description='Compare charges in two sets of molecules and return SMILES strings '
                                                 'of charged centers along with charges in two sets of molecules.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SMILES file with 4 columns: SMILES (for reference), molecule name, the first set '
                             'of protonated SMILES and the second set of protonated SMILES. The file should have '
                             'a header.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output text file.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    pool = Pool(max(min(cpu_count(), args.ncpu), 1))

    with open(args.input) as f:
        with open(args.output, 'wt') as fout:
            header = f.readline().strip().split()
            fout.write(f'smi\tid\tpattern\t{header[2]}_charge\t{header[3]}_charge\n')
            for i, res in enumerate(pool.imap(process_line, f, chunksize=10), 1):
                if res:
                    for item in res:
                        fout.write(item)
                if i % 1000 == 0:
                    sys.stderr.write(f'\rProcessed {i} lines')

            # for line in f:
            #     line = line.strip().split()
            #     res = compare_smiles(line[2], line[3])
            #     for items in res:
            #         fout.write(line[0] + '\t' + line[1] + '\t' + '\t'.join(map(str, items)) + '\n')


if __name__ == '__main__':
    main()

