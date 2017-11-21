#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from read_input import read_input
from multiprocessing import Pool


def get_matrix_best(m, diag=None):
    cids = [c.GetId() for c in m.GetConformers()]
    b = np.full((len(cids), len(cids)), diag, dtype=float)
    for i in range(len(cids)):
        for j in range(i + 1, len(cids)):
            b[i, j] = b[j, i] = AllChem.GetBestRMS(m, m, cids[i], cids[j])
    return b


def filter_conf(data):
    # data is a tuple (mol, mol_name)
    mol, mol_name = data
    if noH:
        mol = Chem.RemoveHs(mol)
    b = get_matrix_best(mol)
    ids = [c.GetId() for c in mol.GetConformers()]
    removed_ids = []
    while b.size and np.nanmin(b) < rms:
        i = np.where(b == np.nanmin(b))[0][0]
        removed_ids.append(ids.pop(i))
        b = np.delete(np.delete(b, i, 0), i, 1)
    for cid in removed_ids:
        mol.RemoveConformer(cid)
    return mol, mol_name


def pool_init(value1, value2):
    global rms
    global noH
    rms = value1
    noH = value2


def main(in_fname, out_fname, rms, noH, ncpu, verbose):
    p = Pool(ncpu, initializer=pool_init, initargs=[rms, noH])
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
                        help='RMS value to filter out conformers. Default: 1.')
    parser.add_argument('-x', '--noH', action='store_true', default=False,
                        help='if set only heavy atoms will be used for RMS calculation.')
    parser.add_argument('-c', '--ncpu', metavar='cpu_number', default=1,
                        help='number of cpus to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "output": out_fname = v
        if o == "rms": rms = float(v)
        if o == "noH": noH = v
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    main(in_fname, out_fname, rms, noH, ncpu, verbose)
