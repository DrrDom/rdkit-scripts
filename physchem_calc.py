#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import sys
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from multiprocessing import Pool, cpu_count


def calc(smi, name):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        hba = rdMolDescriptors.CalcNumHBA(m)
        hbd = rdMolDescriptors.CalcNumHBD(m)
        nrings = rdMolDescriptors.CalcNumRings(m)
        rtb = rdMolDescriptors.CalcNumRotatableBonds(m)
        psa = rdMolDescriptors.CalcTPSA(m)
        logp, mr = rdMolDescriptors.CalcCrippenDescriptors(m)
        mw = rdMolDescriptors._CalcMolWt(m)
        csp3 = rdMolDescriptors.CalcFractionCSP3(m)
        fmf = GetScaffoldForMol(m).GetNumAtoms(onlyHeavy=True) / m.GetNumAtoms(onlyHeavy=True)
        return name, hba, hbd, hba + hbd, nrings, rtb, round(psa, 2), round(logp, 2), round(mr, 2), round(mw, 2), \
               round(csp3, 3), round(fmf, 3)
    else:
        sys.stderr.write('smiles %s cannot be parsed (%s)' % (smi, name))
        return None


def calc_mp(items):
    return calc(*items)


def read_smi(fname, sep="\t"):
    with open(fname) as f:
        for line in f:
            items = line.strip().split(sep)
            if len(items) == 1:
                yield items[0], items[0]
            else:
                yield items[0], items[1]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate some physicochemical parameters with RDKit.')
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input SMILES file. Should contain mol title as a second field.'
                             'Fields are tab-separated. No header.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with calculated physicochemical properties. '
                             'Molecules causing errors will be reported to stderr.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1,
                        help='Number of CPU cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    p = Pool(min(ncpu, cpu_count()))

    with open(out_fname, 'wt') as f:
        f.write('\t'.join(['Name', 'HBA', 'HBD', 'complexity', 'NumRings', 'RTB', 'TPSA', 'logP', 'MR', 'MW']) + '\n')
        for i, res in enumerate(p.imap(calc_mp, read_smi(in_fname), chunksize=100)):
            f.write('\t'.join(map(str, res)) + '\n')
            if verbose and i % 100 == 0:
                sys.stderr.write('\r%i molecules passed' % (i + 1))
                sys.stderr.flush()
