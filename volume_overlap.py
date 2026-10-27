#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from read_input import read_input
from rdkit.Chem import rdShapeHelpers


def read_mols(fname):
    # mol files are read as single-record SDF
    input_format = 'sdf' if fname.lower().endswith('.mol') else None
    for mol, mol_name in read_input(fname, input_format=input_format):
        if mol is None:
            sys.stderr.write(f'Molecule {mol_name} from {fname} cannot be read and was skipped\n')
            continue
        yield mol, mol_name


def calc(ref_fname, query_fnames, output_fname, grid_spacing, ignore_hs, verbose):
    refs = list(read_mols(ref_fname))
    if not refs:
        sys.stderr.write(f'No reference molecules were read from {ref_fname}\n')
        sys.exit(1)

    f = open(output_fname, 'wt') if output_fname is not None else sys.stdout
    i = 0
    try:
        f.write('ref_name\tquery_name\toverlap\n')
        for query_fname in query_fnames:
            for query, query_name in read_mols(query_fname):
                i += 1
                for ref, ref_name in refs:
                    # fraction of the reference volume covered by the query, molecules are not aligned
                    overlap = 1.0 - rdShapeHelpers.ShapeProtrudeDist(ref, query,
                                                                     gridSpacing=grid_spacing,
                                                                     ignoreHs=ignore_hs,
                                                                     allowReordering=False)
                    f.write(f'{ref_name}\t{query_name}\t{round(overlap, 4)}\n')
                if verbose and i % 100 == 0:
                    sys.stderr.write(f'\r{i} query molecules were processed')
    finally:
        if output_fname is not None:
            f.close()
    if verbose:
        sys.stderr.write(f'\r{i} query molecules were processed\n')


def main():
    parser = argparse.ArgumentParser(description='Calculate volume overlap between reference and query molecules '
                                                 'using their input coordinates (no alignment). The overlap is '
                                                 'the fraction of the reference volume covered by a query molecule '
                                                 '(1 - ShapeProtrudeDist(ref, query)).')
    parser.add_argument('-r', '--ref', metavar='FILENAME', required=True, type=str,
                        help='reference molecule(s) in MOL or SDF format. If SDF contains multiple molecules '
                             'each of them will be used as a reference.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str, nargs='+',
                        help='one or more MOL/SDF files with query molecules.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None, type=str,
                        help='output tab-separated text file with columns: reference name, query name, overlap. '
                             'If omitted output will be redirected to STDOUT.')
    parser.add_argument('-g', '--grid_spacing', metavar='NUMERIC', default=0.25, type=float,
                        help='grid spacing used for volume calculation. Default: 0.25.')
    parser.add_argument('--keep_hs', action='store_true', default=False,
                        help='take hydrogens into account. By default hydrogens are ignored.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    calc(args.ref, args.input, args.output, args.grid_spacing, not args.keep_hs, args.verbose)


if __name__ == '__main__':
    main()
