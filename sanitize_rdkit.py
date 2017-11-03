#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import sys
import argparse
from rdkit import Chem
from read_input import read_input


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


def main(in_fname, out_fname, error_fname, input_format, output_format):

    if out_fname is None:
        out_fname = '/dev/stdout'
    if error_fname is None:
        error_fname = '/dev/stderr'

    if output_format == 'sdf':
        wout = Chem.SDWriter(out_fname)
    else:
        print("Error output format.")
        exit()

    werr = Chem.SmilesWriter(error_fname, delimiter='\t', isomericSmiles=True)
    werr.SetProps(['Sanitization_error'])

    for mol, mol_name in read_input(in_fname, input_format=input_format, sanitize=False):

        if mol is not None:
            
            err = Chem.SanitizeMol(mol, catchErrors=True)
            
            if err:

                sys.stderr.write('Error %i sanitizing molecule %s\n' % (err, mol_name))
                sys.stderr.flush()
                mol.SetProp('Sanitization_error', str(err))
                werr.write(mol)

            else:

                try:
                    Chem.Kekulize(mol)
                    double = len(get_unspec_double_bonds(mol))
                    if double:
                        mol.SetProp("Double bonds", str(double))
                    tetra = sum(i[1] == '?' for i in Chem.FindMolChiralCenters(mol, includeUnassigned=True))
                    if tetra:
                        mol.SetProp("Undefined stereocenters", str(tetra))
                    chg = Chem.GetFormalCharge(mol)
                    if chg:
                        mol.SetProp("Charge", str(chg))
                    wout.write(mol)
                except:
                    mol.SetProp('Sanitization_error', 'Kekulization failed')
                    werr.write(mol)

    werr.close()
    wout.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Removes Multi-component compounds and compounds with '
                                                 'non-organic atoms.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SDF format.')
    parser.add_argument('-f', '--input_format', metavar='sdf', required=False, default='sdf',
                        help='input file format. Default: sdf.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=False, default=None,
                        help='output file in SDF format. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-g', '--output_format', metavar='sdf', required=False, default='sdf',
                        help='output file format. Default: sdf.')
    parser.add_argument('-e', '--error_output', metavar='output.smi', required=False, default=None,
                        help='output SMILES file for compounds having sanitization errors. '
                             'If omitted output will be redirected to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "error_output": error_fname = v
        if o == "input_format": input_format = v
        if o == "output_format": output_format = v

    if in_fname == "/dev/stdin":
        in_fname = None

    main(in_fname, out_fname, error_fname, input_format, output_format)


