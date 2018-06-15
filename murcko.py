#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import sys
from multiprocessing import Pool
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmilesFromSmiles


def process_line(line):
    tmp = line.rstrip().split()
    try:
        scaff = MurckoScaffoldSmilesFromSmiles(tmp[0], includeChirality=False)
    except ValueError:
        scaff = ""
    if len(tmp) > 1:
        return tmp[0], tmp[1], scaff
    else:
        return tmp[0], scaff


if __name__ == '__main__':

    pool = Pool(processes=30)

    header = True

    with open(sys.argv[1]) as f:

        if header:
            f.readline()

        print('SMILES\tName\tscaffold')
        for res in pool.imap(process_line, f, chunksize=100):
            print('\t'.join(map(str, res)))

    pool.close()


