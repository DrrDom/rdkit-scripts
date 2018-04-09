#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import numpy as np


def read_thresholds(fname):
    """
    Input file has lines with "prop_name\t1:2:3..."
    :param fname:
    :return:
    """
    d = {}
    with open(fname) as f:
        for line in f:
            items = line.strip().split('\t')
            d[items[0]] = sorted(list(map(float, items[1].split(':'))))
    return d


def load_x(fname):
    x = []
    mol_names = []
    with open(fname) as f:
        var_names = f.readline().strip().split()[1:]
        for line in f:
            items = line.strip().split()
            x.append(list(map(float, items[1:])))
            mol_names.append(items[0])
    return np.array(x), np.array(var_names), np.array(mol_names)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Binning of input values according to specified thresholds.')
    parser.add_argument('-i', '--in', metavar='input.txt', required=True,
                        help='input text file (tab-separated). The first column contains compound names. '
                             'Header contains property names.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with binned values.')
    parser.add_argument('-t', '--thresholds', metavar='thresholds.txt', required=True,
                        help='text file with threshold values.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "thresholds": tr_fname = v

    t = read_thresholds(tr_fname)
    x, var_names, mol_names = load_x(in_fname)

    # keep only vars present in threshold file
    ids = np.in1d(var_names, list(t.keys()))
    x = x[:, ids]
    var_names = var_names[ids]

    for i, v in enumerate(var_names):
        x[:, i] = np.digitize(x[:, i], t[v])

    x = x.astype(int)

    with open(out_fname, 'wt') as f:
        f.write('Name\t' + '\t'.join(var_names) + '\n')
        for i, row in enumerate(x):
            f.write(mol_names[i] + '\t' + '\t'.join(map(str, row)) + '\n')
