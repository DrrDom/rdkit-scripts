#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import re
import sys
from collections import Counter

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition

from read_input import read_input


def read_mols(fname, input_format):
    mols = []
    for mol, mol_name in read_input(fname, input_format=input_format):
        if mol is None:
            sys.stderr.write(f'Molecule {mol_name} cannot be read and was skipped\n')
            continue
        mol = Chem.RemoveHs(mol)
        mol.SetProp('_Name', mol_name)
        mols.append((mol, mol_name))
    return mols


def parse_core(core_str):
    """
    Parse a user-supplied core. Strings with SMARTS-specific syntax ([#6], ;, &, !, ~, comma, $) are parsed
    as SMARTS, all others as SMILES with a fallback to SMARTS. SMILES is preferred because RDKit's SMILES
    parser silently accepts some SMARTS syntax and produces a wrong core, while kekulized SMILES cannot be
    parsed as SMARTS. Dummy atoms with atom map numbers ([*:1]) are treated as attachment points by RDKit.
    """
    if re.search(r'[;&!~,$]|\[#', core_str):
        core = Chem.MolFromSmarts(core_str)
    else:
        core = Chem.MolFromSmiles(core_str)
        if core is None:
            core = Chem.MolFromSmarts(core_str)
    if core is None:
        raise ValueError(f'Core "{core_str}" cannot be parsed as SMILES or SMARTS')
    return core


def find_core(mols, timeout, verbose):
    """
    Identify the common core as the maximum common substructure of all molecules.
    Only complete rings are allowed in the core, ring atoms match ring atoms only and bond orders
    must match exactly (aromatic rings are not matched to aliphatic ones), so the core is a proper scaffold.
    """
    if len(mols) == 1:
        # MCS requires at least two molecules, a single molecule is its own core
        if verbose:
            sys.stderr.write('Only one molecule was read, it is used as the core\n')
        return Chem.Mol(mols[0])
    params = rdFMCS.MCSParameters()
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrderExact
    params.AtomCompareParameters.MatchValences = False
    params.AtomCompareParameters.MatchChiralTag = False
    params.AtomCompareParameters.RingMatchesRingOnly = True
    params.AtomCompareParameters.CompleteRingsOnly = True   # no dangling single ring atoms in the core
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    params.Timeout = timeout
    res = rdFMCS.FindMCS(mols, params)
    if res.canceled:
        sys.stderr.write(f'WARNING: MCS search was stopped by timeout ({timeout} s), '
                         f'the identified core may be not the largest one. '
                         f'Increase --timeout or supply the core explicitly with --core.\n')
    if res.numAtoms == 0:
        raise ValueError('No common core was found for input molecules')
    if verbose:
        sys.stderr.write(f'Core (MCS): {res.smartsString} ({res.numAtoms} atoms, {res.numBonds} bonds)\n')
    return Chem.MolFromSmarts(res.smartsString)


def strip_attachment_points(core):
    """
    Return a copy of the core without terminal dummy atoms (R-group attachment points).
    Such a skeleton can be matched to molecules directly to check how many ways the core fits.
    """
    rw = Chem.RWMol(core)
    to_remove = []
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetDegree() == 1:
            if not a.HasQuery() or a.DescribeQuery().strip().startswith('AtomNull'):
                to_remove.append(a.GetIdx())
    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)
    skeleton = rw.GetMol()
    try:
        skeleton.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(skeleton)
    except Exception:
        pass
    return skeleton


def count_alternative_decompositions(mol, skeleton, max_matches):
    """
    Enumerate all substructure matches of the core skeleton (including symmetry-equivalent ones)
    and count how many of them lead to different assignments of substituents to core positions.
    A value greater than 1 means that R-group positions for this molecule are ambiguous.
    """
    matches = mol.GetSubstructMatches(skeleton, uniquify=False, useChirality=False, maxMatches=max_matches)
    signatures = set()
    for match in matches:
        side_chains = Chem.ReplaceCore(mol, skeleton, match, labelByIndex=True)
        if side_chains is None:
            continue
        # dummy atoms are labelled with the index of a core atom (isotope), shift by one to distinguish index 0
        for a in side_chains.GetAtoms():
            if a.GetAtomicNum() == 0:
                a.SetIsotope(a.GetIsotope() + 1)
        signatures.add(Chem.MolToSmiles(side_chains))
    return len(signatures), len(matches)


def check_ambiguity(mols, core, max_matches):
    skeleton = strip_attachment_points(core)
    ambiguous = []
    for mol, mol_name in mols:
        n_alt, n_matches = count_alternative_decompositions(mol, skeleton, max_matches)
        if n_alt > 1:
            ambiguous.append((mol_name, n_alt))
        if n_matches >= max_matches:
            sys.stderr.write(f'WARNING: {mol_name}: the number of core matches reached the limit '
                             f'({max_matches}), the ambiguity check may be incomplete\n')
    return ambiguous


def decompose(mols, core):
    params = rdRGroupDecomposition.RGroupDecompositionParameters()
    params.removeHydrogensPostMatch = True
    params.removeAllHydrogenRGroups = True   # drop R-groups which are hydrogens for all molecules
    rgd = rdRGroupDecomposition.RGroupDecomposition(core, params)

    matched = []   # indices of matched molecules in the order of rows returned by RDKit
    unmatched = []   # names of molecules which do not match the core
    for i, (mol, mol_name) in enumerate(mols):
        try:
            idx = rgd.Add(mol)
        except Exception as e:
            sys.stderr.write(f'WARNING: {mol_name} cannot be decomposed: {e}\n')
            idx = -1
        if idx < 0:
            unmatched.append(mol_name)
        else:
            matched.append(i)

    columns = {}
    if matched:
        rgd.Process()
        columns = rgd.GetRGroupsAsColumns(asSmiles=True)
    return matched, unmatched, columns


def canonical_core(smi):
    """
    Canonical core SMILES without stereo (the core is realized from a particular molecule by RDKit,
    its stereo is not a property of the common core).
    """
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def write_table(output_fname, mols, matched, columns):
    rgroup_names = sorted((k for k in columns if k != 'Core'), key=lambda x: int(x[1:]))

    cores = Counter(canonical_core(smi) for smi in columns.get('Core', []))
    if len(cores) > 1:
        sys.stderr.write('WARNING: core representations differ between molecules '
                         '(e.g. due to different charges or hydrogens on core atoms), '
                         'the most common one is written in the header. All found cores: ' +
                         ', '.join(cores.keys()) + '\n')
    core_smi = cores.most_common(1)[0][0] if cores else ''

    rows = {mol_idx: [columns[r][row_idx] for r in rgroup_names] for row_idx, mol_idx in enumerate(matched)}

    f = open(output_fname, 'wt') if output_fname is not None else sys.stdout
    try:
        f.write(f'# core: {core_smi}\n')
        f.write('\t'.join(['name'] + rgroup_names) + '\n')
        for i, (mol, mol_name) in enumerate(mols):
            values = rows.get(i, [''] * len(rgroup_names))
            f.write('\t'.join([mol_name] + values) + '\n')
    finally:
        if output_fname is not None:
            f.close()


def main():
    parser = argparse.ArgumentParser(description='R-group decomposition of molecules sharing a common core. '
                                                 'The core is identified automatically as the maximum common '
                                                 'substructure (complete rings only) of all input molecules, '
                                                 'unless it is supplied explicitly. The output is a tab-separated '
                                                 'table: molecule name and R-groups (R1, R2, ...) as SMILES with '
                                                 'attachment points labelled [*:N]; unsubstituted positions are '
                                                 'reported as hydrogens ([H][*:N]). The first line of the output '
                                                 'starts with # and contains the core SMILES with labelled '
                                                 'R-group positions. If the core can be matched to a molecule in '
                                                 'several non-equivalent ways (symmetric cores, e.g. naphthalene) '
                                                 'a warning is printed, because the alignment of R-groups between '
                                                 'such molecules is ambiguous and is chosen by RDKit heuristically. '
                                                 'Note: substituents common to all molecules become part of the '
                                                 'automatically detected core.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None, type=str,
                        help='input SMILES (tab-separated, names in the second column) or SDF file '
                             '(molecule titles are used as names). If omitted SMILES will be read from STDIN.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None, type=str,
                        help='output tab-separated text file. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-c', '--core', metavar='SMILES/SMARTS', required=False, default=None, type=str,
                        help='core to use instead of the automatically detected one. Attachment points may be '
                             'labelled with atom maps, e.g. [*:1]c1ccc([*:2])cc1; otherwise R-groups are detected '
                             'at all core positions. Strings with SMARTS-specific syntax ([#6], ;, &, !, ~, comma) '
                             'are parsed as SMARTS, all others as SMILES.')
    parser.add_argument('-t', '--timeout', metavar='INTEGER', required=False, default=60, type=int,
                        help='timeout (seconds) of the maximum common substructure search.')
    parser.add_argument('--max_matches', metavar='INTEGER', required=False, default=10000, type=int,
                        help='maximum number of core matches per molecule enumerated for the ambiguity check.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print additional information to STDERR.')
    args = parser.parse_args()

    input_fname = None if args.input == '/dev/stdin' else args.input
    input_format = 'smi' if input_fname is None else None
    mols = read_mols(input_fname, input_format)
    if not mols:
        sys.stderr.write(f'No molecules were read from {args.input if input_fname else "STDIN"}\n')
        sys.exit(1)
    if args.verbose:
        sys.stderr.write(f'{len(mols)} molecules were read\n')

    names = Counter(name for _, name in mols)
    duplicated = [name for name, n in names.items() if n > 1]
    if duplicated:
        sys.stderr.write('WARNING: duplicated molecule names, the table will contain rows with the same name: ' +
                         ', '.join(duplicated) + '\n')

    try:
        if args.core is not None:
            core = parse_core(args.core)
        else:
            core = find_core([m for m, _ in mols], args.timeout, args.verbose)
    except ValueError as e:
        sys.stderr.write(f'{e}\n')
        sys.exit(1)

    ambiguous = check_ambiguity(mols, core, args.max_matches)
    if ambiguous:
        sys.stderr.write(f'WARNING: the core can be matched in several non-equivalent ways to {len(ambiguous)} '
                         f'of {len(mols)} molecules, so the assignment of their substituents to R-group positions '
                         f'is ambiguous (e.g. for symmetric cores like naphthalene). RDKit aligned them '
                         f'heuristically by minimising the number of distinct R-groups; check the result or '
                         f'supply a labelled core with --core. Affected molecules (number of alternative '
                         f'assignments): ' + ', '.join(f'{name} ({n})' for name, n in ambiguous) + '\n')

    matched, unmatched, columns = decompose(mols, core)
    if unmatched:
        sys.stderr.write(f'WARNING: {len(unmatched)} molecules do not match the core and have empty R-groups '
                         f'in the table: ' + ', '.join(unmatched) + '\n')
    if matched and not [k for k in columns if k != 'Core']:
        sys.stderr.write('WARNING: no R-groups were found, all matched molecules are identical to the core\n')
    if args.verbose:
        sys.stderr.write(f'{len(matched)} molecules were decomposed, '
                         f'{len([k for k in columns if k != "Core"])} R-group positions were found\n')

    write_table(args.output, mols, matched, columns)


if __name__ == '__main__':
    main()
