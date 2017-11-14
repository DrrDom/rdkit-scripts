#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
import pickle
import numpy as np
from rdkit.Chem import AllChem_2   # works only in rdk_2017_3 env!!
from read_input import read_input
from multiprocessing import Pool


def filter_conf(data):
    # data is a tuple (mol, mol_name)
    mol, mol_name = data
    ids = [c.GetId() for c in mol.GetConformers()]
    a = AllChem_2.GetConformerRMSMatrix(mol)
    b = np.full((len(ids), len(ids)), None, dtype=float)
    b[np.tril_indices(len(ids), -1)] = a
    removed_ids = []
    while b.size and np.nanmin(b) < rms:
        i = np.where(b == np.nanmin(b))[0][0]
        removed_ids.append(ids.pop(i))
        b = np.delete(np.delete(b, i, 0), i, 1)
    for cid in removed_ids:
        mol.RemoveConformer(cid)
    return mol, mol_name


def pool_init(value):
    global rms
    rms = value


def main(in_fname, out_fname, rms, ncpu, verbose):
    p = Pool(ncpu, initializer=pool_init, initargs=[rms])
    with open(out_fname, 'wb') as fout:
        for i, res in enumerate(p.imap_unordered(filter_conf, read_input(in_fname), chunksize=20), 1):
            pickle.dump(res, fout, -1)
            if verbose and (i + 1) % 1000 == 0:
                sys.stderr.write('\rCompounds processed: %i' % (i + 1))
    sys.stderr.write('\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter conformers by RMS.')
    parser.add_argument('-i', '--input', metavar='input.pkl', required=True,
                        help='input file in PKL format containing pickled molecules with multiple conformations.')
    parser.add_argument('-o', '--output', metavar='output.pkl', required=True,
                        help='output file in PKL format with conformers filtered out according to specified RMS value.')
    parser.add_argument('-r', '--rms', required=False, default=1,
                        help='RMS value to filter out conformers.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "rms": rms = float(v)
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    main(in_fname, out_fname, rms, ncpu, verbose)
