#!/usr/bin/env python3

__author__ = 'pavel'

import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import product
from copy import deepcopy
from multiprocessing import Pool, cpu_count
from read_input import read_input


def prep_input(fname, id_field_name, tetrahedral, double_bond, max_undef):
    for mol, mol_name in read_input(fname, id_field_name=id_field_name):
        yield mol, mol_name, tetrahedral, double_bond, max_undef


def map_enumerate_stereo(args):
    return enumerate_stereo(*args)


def get_unspec_double_bonds(m):

    def check_nei_bonds(bond):
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        a1_bonds_single = [b.GetBondType() == Chem.BondType.SINGLE for b in a1.GetBonds() if b.GetIdx() != bond.GetIdx()]
        a2_bonds_single = [b.GetBondType() == Chem.BondType.SINGLE for b in a2.GetBonds() if b.GetIdx() != bond.GetIdx()]

        # if there are two identical substituents in one side then the bond is unsteric (no stereoisomers possible)
        ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
        a1_nei = [a.GetIdx() for a in a1.GetNeighbors() if a.GetIdx() != a2.GetIdx()]
        if len(a1_nei) == 2 and \
                all(m.GetBondBetweenAtoms(i, a1.GetIdx()).GetBondType() == Chem.BondType.SINGLE for i in a1_nei) and \
                ranks[a1_nei[0]] == ranks[a1_nei[1]]:
            return False
        a2_nei = [a.GetIdx() for a in a2.GetNeighbors() if a.GetIdx() != a1.GetIdx()]
        if len(a2_nei) == 2 and \
                all(m.GetBondBetweenAtoms(i, a2.GetIdx()).GetBondType() == Chem.BondType.SINGLE for i in a2_nei) and \
                ranks[a2_nei[0]] == ranks[a2_nei[1]]:
            return False

        # if list is empty this is a terminal atom, e.g. O in C=O
        if a1_bonds_single and a2_bonds_single and \
                all(a1_bonds_single) and all(a2_bonds_single):
            return True
        else:
            return False

    res = []
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and \
           b.GetStereo() == Chem.BondStereo.STEREONONE and \
           (not b.IsInRing() or not (b.IsInRingSize(3) or b.IsInRingSize(4) or b.IsInRingSize(5) or b.IsInRingSize(6) or b.IsInRingSize(7))) and \
           check_nei_bonds(b):
            res.append(b.GetIdx())
    return res


def set_double_bond_stereo(bond, value):
    # value can be 1 or 0

    def set_bond(fixed_bonds, changed_bonds, bond_value):

        if bond_value:
            if any(b.GetBondDir() == Chem.rdchem.BondDir.ENDDOWNRIGHT for b in fixed_bonds):
                changed_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)
            elif any(b.GetBondDir() == Chem.rdchem.BondDir.ENDUPRIGHT for b in fixed_bonds):
                changed_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)
        else:
            if any(b.GetBondDir() == Chem.rdchem.BondDir.ENDDOWNRIGHT for b in fixed_bonds):
                changed_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)
            elif any(b.GetBondDir() == Chem.rdchem.BondDir.ENDUPRIGHT for b in fixed_bonds):
                changed_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)

    a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
    a1_bonds = tuple(b for b in a1.GetBonds() if b.GetIdx() != bond.GetIdx())
    a2_bonds = tuple(b for b in a2.GetBonds() if b.GetIdx() != bond.GetIdx())

    if a1_bonds and a2_bonds:   # if double bond of carbonyl one atom does not have any other bonds - skip

        if all(b.GetBondDir() == Chem.BondDir.NONE for b in a1_bonds) and \
           all(b.GetBondDir() == Chem.BondDir.NONE for b in a2_bonds):

            a1_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)
            if value:
                a2_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)
            else:
                a2_bonds[0].SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)

        elif all(b.GetBondDir() == Chem.BondDir.NONE for b in a1_bonds):
            set_bond(a2_bonds, a1_bonds, value)
        else:
            set_bond(a1_bonds, a2_bonds, value)


def enumerate_double_bond_stereo(mol):
    bonds = get_unspec_double_bonds(mol)
    res = []
    for p in product([0, 1], repeat=len(bonds)):
        m = deepcopy(mol)
        for value, bond in zip(p, bonds):
            set_double_bond_stereo(m.GetBondWithIdx(bond), value)
        Chem.AssignStereochemistry(m, force=True, cleanIt=True)
        res.append(m)
    return res


def enumerate_tetrahedral_stereo(mol):

    output = []
    chiral_undef = tuple(i[0] for i in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if i[1] == '?')

    # without Hs sometimes wrong stereochemistry can be optimized without errors CC(C)=CCN1CCC2(C)c3cc(O)ccc3CC1C2C
    mol = Chem.AddHs(mol)

    num = 0
    for p in product([Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW], repeat=len(chiral_undef)):
        for tag, i in zip(p, chiral_undef):
            mol.GetAtomWithIdx(i).SetChiralTag(tag)
        AllChem.EmbedMolecule(mol)
        try:
            if AllChem.UFFHasAllMoleculeParams(mol):
                AllChem.UFFOptimizeMolecule(mol, maxIters=10)
                num += 1
                output.append(deepcopy(Chem.RemoveHs(mol)))
            else:
                sys.stderr.write('No UFF parameters for %s\n' % (Chem.MolToSmiles(Chem.RemoveHs(mol))))
        except ValueError:
            continue

    return output


def enumerate_stereo(mol, mol_name, tetrahedral, double_bond, max_undef):

    if max_undef != -1:
        undef = 0
        if tetrahedral:
            undef += sum(i[1] == '?' for i in Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        if double_bond:
            undef += len(get_unspec_double_bonds(mol))
        if undef > max_undef:
            return []

    mols = [mol]
    if double_bond:
        mols = enumerate_double_bond_stereo(mol)
    if tetrahedral:
        tmp = []
        for m in mols:
            tmp.extend(enumerate_tetrahedral_stereo(m))
        mols = tmp

    # filter possible duplicates and keep the order, not sure whether it is really needed
    output = []
    for m in mols:
        smi = Chem.MolToSmiles(m, isomericSmiles=True)
        if smi not in output:
            output.append(smi)

    return [(smi, "%s_%i" % (mol_name, i+1)) for i, smi in enumerate(output)]


def main_params(in_fname, out_fname, tetrahedral, double_bond, max_undef, id_field_name, ncpu, verbose):

    if out_fname is not None:
        fout = open(out_fname, 'wt')

    nprocess = min(cpu_count(), max(ncpu, 1))
    p = Pool(nprocess)

    try:
        for i, res in enumerate(p.imap_unordered(map_enumerate_stereo,
                                                 prep_input(in_fname, id_field_name, tetrahedral, double_bond, max_undef),
                                                 chunksize=10)):
            if out_fname is None:
                for smi, mol_name in res:
                    print(smi + '\t' + mol_name)
                sys.stdout.flush()
            else:
                for smi, mol_name in res:
                    fout.write('%s\t%s\n' % (smi, mol_name))
            if verbose and i % 100 == 0:
                sys.stderr.write('\r%i molecules were processed' % i)
                sys.stderr.flush()
    finally:
        p.close()

    if out_fname is not None:
        fout.close()


def main():

    parser = argparse.ArgumentParser(description='Generation of stereoisomers with RDKit.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=True,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional.')
    parser.add_argument('-o', '--output', metavar='output.smi', required=False, default=None,
                        help='output file in SMILES format. If omitted output will be made in STDOUT.')
    parser.add_argument('-t', '--tetrahedral', required=False, action='store_true', default=False,
                        help='generate stereoisomers for unspecified tetrahedral centers.')
    parser.add_argument('-d', '--double_bond', required=False, action='store_true', default=False,
                        help='generate stereoisomers for unspecified double bonds.')
    parser.add_argument('-u', '--max_undef', metavar='INTEGER', default=-1,
                        help='maximum allowed number of unspecified stereocenters and/or double bonds. '
                             'if compound contains greater number of them it will be discarded. '
                             'Default: all possible stereoisomers will be enumerated '
                             '(beware of combinatorial explosion).')
    parser.add_argument('-f', '--id_field_name', metavar='field_name', default=None,
                        help='SDF input: field name of compound ID. '
                             'If omitted molecule titles will be used or SMILES string as name.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "id_field_name": id_field_name = v
        if o == "ncpu": ncpu = int(v)
        if o == "max_undef": max_undef = int(v)
        if o == "tetrahedral": tetrahedral = v
        if o == "double_bond": double_bond = v
        if o == "verbose": verbose = v

    if not tetrahedral and not double_bond:
        print("You should specify at least one option -t or -d. Revise you command line arguments.")
        exit()

    main_params(in_fname, out_fname, tetrahedral, double_bond, max_undef, id_field_name, ncpu, verbose)


if __name__ == '__main__':
    main()
