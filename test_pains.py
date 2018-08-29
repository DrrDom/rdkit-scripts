#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import sys
import argparse
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from multiprocessing import Pool, cpu_count


def read_smiles(fname, smiles_col, names_col, header):
    with open(fname) as f:
        if header:
            f.readline()
        for line in f:
            tmp = line.split()
            yield tmp[smiles_col], tmp[names_col]


def check_smiles(input):
    smi, name = input
    res = []
    try:
        mol = Chem.MolFromSmiles(smi)
        for entry in catalog.GetFilterMatches(mol):
            res.append("\t".join(map(str, (smi, name, entry.filterMatch))))
    except:
        res.append(smi + '\t' + name + '\t' + 'parse error')
    return res


def init():
    global catalog
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog(params)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test input molecules for PAINS.')
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input SMILES file. Fields are whitespace-separated.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with SMILES, names and found PAINS.')
    parser.add_argument('-s', '--smiles_col', default=0,
                        help='column number in the input file which contains SMILES. Default: 0.')
    parser.add_argument('-n', '--names_col', default=1,
                        help='column number in the input file which contains molecule titles. Default: 1.')
    parser.add_argument('--header', action='store_true', default=False,
                        help='if specified the first line (header) of the input will be omitted.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1,
                        help='number of CPU cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "smiles_col": smiles_col = int(v)
        if o == "names_col": names_col = int(v)
        if o == "header": header = v
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    p = Pool(min(ncpu, cpu_count()), initializer=init)

    with open(out_fname, 'wt') as f:
        for i, res in enumerate(p.imap_unordered(check_smiles, read_smiles(in_fname, smiles_col, names_col, header), chunksize=100)):
            if res:
                f.write('\n'.join(res) + '\n')
            if verbose and i % 10000 == 0:
                sys.stderr.write('\r%i molecules passed' % (i + 1))
                sys.stderr.flush()

