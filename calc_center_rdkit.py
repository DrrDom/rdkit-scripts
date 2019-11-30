#!/usr/bin/env python2

__author__ = 'Pavel Polishchuk'


from rdkit import Chem
import argparse
import sys


def read_pdbqt(fname, sanitize, removeHs):
    mols = []
    with open(fname) as f:
        pdb_block = f.read().split('MODEL ')
        for j, block in enumerate(pdb_block[1:]):
            m = Chem.MolFromPDBBlock('\n'.join([i[:66] for i in block.split('\n')]),
                                     sanitize=sanitize,
                                     removeHs=removeHs)
            if m is None:
                sys.stderr.write(f'The pose #{j+1} cannot be read from {fname}\n')
            else:
                mols.append(m)
    return mols


def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)


def calc_center(coord):
    # coord tuple of tuples ((x1,y1,z1), (x2,y2,z2), ...)
    if len(coord) == 0:
        return None
    min_x, min_y, min_z = coord[0]
    max_x, max_y, max_z = coord[0]
    for item in coord[1:]:
        if item[0] < min_x:
            min_x = item[0]
        if item[0] > max_x:
            max_x = item[0]
        if item[1] < min_y:
            min_y = item[1]
        if item[1] > max_y:
            max_y = item[1]
        if item[2] < min_z:
            min_z = item[2]
        if item[2] > max_z:
            max_z = item[2]
    output = [(max_x + min_x) / 2, (max_y + min_y) / 2, (max_z + min_z) / 2]
    output = (round(i, 3) for i in output)
    return output


def main_params(input_fnames, output_fname):

    if output_fname is not None:
        sys.stdout = open(output_fname, 'wt')

    for fname in input_fnames:
        if fname[-4:].lower() == '.pdb':
            f = Chem.MolFromPDBFile
        elif fname[-5:].lower() == '.mol2':
            f = Chem.MolFromMol2File
        elif fname[-6:].lower() == '.pdbqt':
            f = read_pdbqt
        else:
            continue

        m = f(fname, sanitize=False, removeHs=False)

        if m is not None:

            if isinstance(m, list):  # pdbqt
                for item in m:
                    center = calc_center(get_coord(item))
                    if center is not None:
                        print(fname + '\t' + '\t'.join(map(str, center)))

            else:  # pdb/mol2

                center = calc_center(get_coord(m))
                if center is not None:
                    print(fname + '\t' + '\t'.join(map(str, center)))

        else:

            sys.stderr.write(f'Molecule cannot be read from {fname}\n')


def main():

    parser = argparse.ArgumentParser(description='Calc center of all supplied coordinates in a file.')
    parser.add_argument('-i', '--input', metavar='input_molecule', required=True, nargs='*',
                        help='input files (PDB/PDBQT/MOL2).')
    parser.add_argument('-o', '--output', metavar='output.txt',
                        help='output text file. If omitted output will be in stdout.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": input_fnames = v
        if o == "output": output_fname = v

    main_params(input_fnames, output_fname)


if __name__ == '__main__':
    main()
