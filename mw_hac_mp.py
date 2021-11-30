#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import sys
from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt


def process_line(line):
    tmp = line.strip().split()
    m = Chem.MolFromSmiles(tmp[0])
    if m:
        mw = CalcExactMolWt(m)
        hac = m.GetNumHeavyAtoms()
        return tmp[0], tmp[1], mw, hac
    else:
        return None


def main():

    pool = Pool(processes=30)
    
    header = True
    
    with open(sys.argv[1]) as f:
        
        if header:
            f.readline()

        print('SMILES\tName\tmw\thac')
        for res in pool.imap(process_line, f, chunksize=100):
            if res:
                print('\t'.join(map(str, res)))

    pool.close()


if __name__ == '__main__':
    main()
