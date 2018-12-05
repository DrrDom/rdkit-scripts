#!/usr/bin/env python2

__author__ = 'Pavel Polishchuk'


from rdkit import Chem
import argparse
import sys


def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)


def rmsd(mol, ref):
    match_indeces = mol.GetSubstructMatches(ref, uniquify=False)
    if not match_indeces:
        return None
    else:
        ref_coord = get_coord(ref)
        mol_coord = get_coord(mol, match_indeces[0])
        s = 0
        for r, m in zip(ref_coord, mol_coord):
            s += ((r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2) ** 0.5
        return round(s / len(ref_coord), 3)


def main_params(input_fnames, output_fname, ref_name):
    ref = Chem.MolFromMol2File(ref_name)

    if output_fname is not None:
        sys.stdout = open(output_fname, 'wt')

    for in_fname in input_fnames:
        mol = Chem.MolFromMol2File(in_fname)
        if mol is None:
            print('%s\t%s' % (in_fname, 'Cannot read structure'))
        else:
            mol_rmsd = rmsd(mol, ref)
            if mol_rmsd is not None:
                print('%s\t%s' % (in_fname, str(mol_rmsd)))
            else:
                print('%s\t%s' % (in_fname, 'No matches'))


def main():

    parser = argparse.ArgumentParser(description='Calc RMSD between ref molecules and docked poses.')
    parser.add_argument('-i', '--input', metavar='input.mol2', required=True, nargs='*',
                        help='input mol2 files to compare with a reference molecule.')
    parser.add_argument('-r', '--reference', metavar='ref.mol2', required=True,
                        help='reference molecule (from X-ray complex structure).')
    parser.add_argument('-o', '--output', metavar='output.txt',
                        help='output text file. If omitted output will be in stdout.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fnames = v
        if o == "output": output_fname = v
        if o == "reference": ref_name = v

    main_params(input_fnames, output_fname, ref_name)


if __name__ == '__main__':
    main()
