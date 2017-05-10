#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'


import argparse
from rdkit import Chem
from read_input import read_input


# list excludes transition metals (they can form complices) and Mg, Ca,... - they can form complices
good_elm = {'H', 'C', 'O', 'N', 'P', 'Cl', 'F', 'Br', 'I', 'S',
               'B', 'Si', "Cs", "Li",
               "Na", "K", "Rb",  "Al", "Ga", "In","Ti", "Nh",
               "Ge", "Sn", "Pb", "Fl", "As", "Sb", "Bi", "Se",
               "Te", "At", "Ts"}


def main_params(in_fname, good_fname, bad_fname):

    if good_fname is not None:
        wgood = Chem.SDWriter(good_fname)
    if bad_fname is not None:
        wbad = Chem.SDWriter(bad_fname)

    for mol, mol_name in read_input(in_fname, input_format='sdf', sanitize=False):

        if mol is not None:
            
            errors = False

            n = len(Chem.GetMolFrags(mol))
            # chg = Chem.GetFormalCharge(mol)
            elm = set(a.GetSymbol() for a in mol.GetAtoms())
            elm_diff = elm - good_elm

            if 'C' not in elm:  # no carbon
                mol.SetProp('Eval:no carbons', '1')
                errors = True
            if elm_diff:        # any illegible atom
                mol.SetProp('Eval:illegible atoms', ', '.join(sorted(list(elm_diff))))
                errors = True
            if n > 1:           # multi-component
                if len(set(Chem.MolToSmiles(frag, isomericSmiles=True) for frag in Chem.GetMolFrags(mol, asMols=True))) > 1:  # number of unique components
                    mol.SetProp('Eval:multi-component', str(n))
                    errors = True
                else:           # if all components are identical keep one component as a mol with all fields
                    tmp_mol = Chem.GetMolFrags(mol, asMols=True)[0]
                    for prop_name in mol.GetPropNames(includePrivate=True):
                        tmp_mol.SetProp(prop_name, mol.GetProp(prop_name))
                    mol = tmp_mol

            if errors:
                if bad_fname is not None:
                    wbad.write(mol)
            else:
                wgood.write(mol)


def main():

    parser = argparse.ArgumentParser(description='Removes Multi-component compounds and compounds with '
                                                 'non-organic atoms.')
    parser.add_argument('-i', '--input', metavar='input.sdf', required=False, default=None,
                        help='input file in SDF or SMILES format. SMILES input should have no header, '
                             'the first column is SMILES string and the second column with ID is optional. '
                             'If omitted STDIN will be read as SDF format.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True,
                        help='output file in SDF format.')
    parser.add_argument('-d', '--discarded', metavar='output.smi', required=False, default=None,
                        help='output file for discarded compounds in SDF or SMILES format. '
                             'If omitted the file will not be created and no output will be.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": good_fname = v
        if o == "discarded": bad_fname = v

    if in_fname == "/dev/stdin":
        in_fname = None

    main_params(in_fname, good_fname, bad_fname)


if __name__ == '__main__':
    main()
