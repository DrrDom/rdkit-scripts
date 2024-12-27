#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
from functools import partial
import sys
from read_input import read_input
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from multiprocessing import Pool, cpu_count


def enum_stereo(item, max_isomers, use_embedding, no_suffix):
    mol, mol_name = item
    opts = StereoEnumerationOptions(tryEmbedding=use_embedding, maxIsomers=max_isomers)
    output = []
    # this is a workaround for rdkit issue - if a double bond has STEREOANY it will cause errors at
    # stereoisomer enumeration, we replace STEREOANY with STEREONONE in these cases
    try:
        isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    except RuntimeError:
        for bond in mol.GetBonds():
            if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                bond.SetStereo(Chem.BondStereo.STEREONONE)
        isomers = tuple(EnumerateStereoisomers(mol,options=opts))
    for i, m in enumerate(isomers, 1):
        name = mol_name if no_suffix else f'{mol_name}_{i}'
        output.append((Chem.MolToSmiles(m, isomericSmiles=True), name))
    return output


def enum_stereoisomers(input_fname, output_fname, max_isomers, use_embedding, no_suffix, ncpu, verbose):

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    with open(output_fname, 'wt') as f:
        for i, items in enumerate(pool.imap(partial(enum_stereo,
                                                    max_isomers=max_isomers,
                                                    use_embedding=use_embedding,
                                                    no_suffix=no_suffix),
                                            read_input(input_fname)), 1):
            for smi, mol_name in items:
                f.write(f'{smi}\t{mol_name}\n')
            if verbose and i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were processed')

    if verbose:
        sys.stderr.write(f'\r{i} molecules were processed\n')


def main():
    parser = argparse.ArgumentParser(description='INSERT DESCRIPTION HERE.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF or SMILES file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SMILES file.')
    parser.add_argument('-m', '--max_isomers', metavar='INTEGER', default=1, type=int,
                        help='maximum number of stereoisomers generated. Default: 1.')
    parser.add_argument('-e', '--embed', action='store_true', default=False,
                        help='use embedding to check validity of stereoisomers.')
    parser.add_argument('-x', '--no_suffix', action='store_true', default=False,
                        help='do not append suffix with a sequential number of an isomer to molecule name.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    enum_stereoisomers(args.input, args.output, args.max_isomers, args.embed, args.no_suffix, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
