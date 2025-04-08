#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import os
import re
import sys

from itertools import product

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdMolAlign import GetBestRMS

from read_input import read_pdbqt, read_input


def get_coord(mol, indices=None):
    if indices is None:
        indices = tuple(range(mol.GetNumAtoms()))
    output = []
    for atom_id in indices:
        pos = mol.GetConformer().GetAtomPosition(atom_id)
        output.append((pos.x, pos.y, pos.z))
    return tuple(output)


def rmsd(mol, ref, chirality, align):

    def rmsd_calc(r_coord, m_coord):
        s = 0
        for r, m in zip(r_coord, m_coord):
            s += (r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2
        s = (s / len(r_coord)) ** 0.5
        return s

    match_indices = mol.GetSubstructMatches(ref, uniquify=False, useChirality=chirality, maxMatches=10000)
    min_rmsd = float('inf')
    if not match_indices:
        mcs = rdFMCS.FindMCS([mol, ref], threshold=1.0,
                             ringMatchesRingOnly=False, completeRingsOnly=False,
                             matchChiralTag=chirality)
        if not mcs:
            return None
        patt = Chem.MolFromSmarts(mcs.smartsString)
        refMatch, molMatch = ref.GetSubstructMatches(patt, uniquify=False), \
                             mol.GetSubstructMatches(patt, uniquify=False)
        if not align:
            for ids_ref in refMatch:
                for ids_mol in molMatch:
                    ref_coord = get_coord(ref, ids_ref)
                    mol_coord = get_coord(mol, ids_mol)
                    s = rmsd_calc(ref_coord, mol_coord)
                    if s < min_rmsd:
                        min_rmsd = s
        else:
            map_ids = [list(zip(ids_mol, ids_ref)) for ids_mol, ids_ref in product(molMatch, refMatch)]
            min_rmsd = GetBestRMS(mol, ref, map=map_ids)

    else:
        if not align:
            ref_coord = get_coord(ref)
            for ids in match_indices:
                mol_coord = get_coord(mol, ids)
                s = rmsd_calc(ref_coord, mol_coord)
                if s < min_rmsd:
                    min_rmsd = s
        else:
            ref_ids = list(range(ref.GetNumAtoms()))
            map_ids = [list(zip(ids, ref_ids)) for ids in match_indices]
            min_rmsd = GetBestRMS(mol, ref, map=map_ids)

    return round(min_rmsd, 3)


def main_params(input_fnames, input_smi, output_fname, ref_name, refsmi, chirality, regex, align):

    if ref_name.lower().endswith('.mol2'):
        ref = Chem.MolFromMol2File(ref_name, removeHs=True)
    elif ref_name.lower().endswith('.pdbqt'):
        ref = read_pdbqt(ref_name, refsmi, removeHs=True)[0]
    elif ref_name.lower().endswith('.mol'):
        ref = Chem.MolFromMolFile(ref_name, removeHs=True)
    elif ref_name.lower().endswith('.sdf'):
        ref = {m.GetProp('_Name'): m for m in Chem.SDMolSupplier(ref_name) if m}
    else:
        sys.stderr.write('Wrong format of the reference file. Only MOL2, PDBQT and SDF files are allowed.\n')
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
                    sys.stderr.write(
                        f'Line "{line}" in input smiles does not have two fields - SMILES and mol name. Skipped.\n')

    for in_fname in input_fnames:

        if in_fname.lower().endswith('.mol2'):
            mols = [Chem.MolFromMol2File(in_fname)]
        elif in_fname.lower().endswith('.pdbqt') or in_fname.endswith('.pdbqt_out'):
            if regex is not None:
                mols = read_pdbqt(in_fname, smis[re.search(regex, os.path.basename(in_fname)).group()], removeHs=True)
            else:
                mols = read_pdbqt(in_fname, smis[os.path.splitext(os.path.basename(in_fname))[0]], removeHs=True)
        elif in_fname.lower().endswith('.sdf'):
            mols = [mol for mol, mol_name in read_input(in_fname)]
        else:
            sys.stderr.write(f'Wrong format of the input file - {in_fname}. '
                             f'Only MOL2, PDBQT and SDF files are allowed.\n')
            raise ValueError

        for i, mol in enumerate(mols, 1):
            if mol is None:
                print(f'{in_fname}\t{i}\tCannot read structure')
            else:
                # assign printed mol name from mol object or file name
                mol_name = mol.GetProp('_Name')
                if not mol_name:
                    if regex is not None:
                        mol_name = re.search(regex, os.path.basename(in_fname)).group()
                    else:
                        mol_name = os.path.basename(in_fname)
                # assign ref mol object
                if isinstance(ref, dict):
                    try:
                        refmol = ref[mol.GetProp('_Name')]
                    except KeyError:
                        sys.stderr.write(f'Molecule with name {mol.GetProp("_Name")} is not available '
                                         f'in the reference SDF file\n')
                        print(f'{mol_name}\t{i}\tNo matches')
                        refmol = None
                else:
                    refmol = ref

                if refmol is not None:
                    mol_rmsd = rmsd(mol, refmol, chirality, align)
                    if mol_rmsd is not None:
                        print(f'{mol_name}\t{i}\t{mol_rmsd}')
                    else:
                        print(f'{mol_name}\t{i}\tNo matches')


def main():
    parser = argparse.ArgumentParser(description='''Calc RMSD between a reference molecule and docked poses.
                                                 If reference molecule is not substructure of the docked molecule
                                                 maximum common substructure is used.''')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, nargs='*',
                        help='input MOL2/PDBQT/SDF file(s) to compare with a reference molecule or molecules.')
    parser.add_argument('--input_smi', metavar='FILENAME', required=False, default=None,
                        help='SMILES of input molecules if they are in PDBQT format. No header. Space- or '
                             'tab-separated. Molecule names should correspond to PDBQT file names. '
                             'DEVELOPERS, STOP USING PDB FORMAT FOR DOCKING PROGRAMS!')
    parser.add_argument('-r', '--reference', metavar='FILENAME', required=True,
                        help='reference molecule (from X-ray complex structure) in MOL/MOL2/PDBQT/SDF format. '
                             'Multiple reference molecules can be supplied in SDF file, in this case molecules '
                             'will be matched by their names.')
    parser.add_argument('-s', '--refsmi', metavar='SMILES or FILENAME', required=False, default=None,
                        help='SMILES of the reference molecule. This can be a string (do not forget to escape '
                             'special characters) or a file. This is required only for PDBQT input to '
                             'assign bond orders. DEVELOPERS, STOP USING PDB FORMAT FOR DOCKING PROGRAMS!')
    parser.add_argument('--regex', metavar='REGEX', required=False, default=None,
                        help='Use it if there are complex names of pdbqt files. '
                             'Use regex search to establish a relationship between reference smiles name and pdbqt '
                             'filename. If None filename of pdbqt file will be taken to find reference smiles '
                             'Examples: MOLID[0-9]* or .*(?=_out.pdbqt)')
    parser.add_argument('-o', '--output', metavar='FILENAME',
                        help='output text file. If omitted output will be in stdout.')
    parser.add_argument('-a', '--align', action='store_true', default=False,
                        help='align molecules to calculate minimum RMSD. By default there is no alignment, that is '
                             'suitable for assessing of docking outputs.')
    parser.add_argument('-x', '--nochirality', action='store_true', default=False,
                        help='choose this option if you want to omit matching chirality in substructure search. '
                             'By default chirality is considered.')

    args = parser.parse_args()
    if (args.refsmi is not None) and (args.refsmi.lower().endswith('.smi') or args.refsmi.lower().endswith('.smiles')):
        with open(args.refsmi) as inp:
            refsmi = inp.read().strip()
    else:
        refsmi = args.refsmi

    main_params(args.input, args.input_smi, args.output, args.reference, refsmi, not args.nochirality, args.regex, args.align)


if __name__ == '__main__':
    main()
