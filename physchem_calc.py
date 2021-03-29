#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import sys
import argparse
from argparse import RawTextHelpFormatter
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, QED, FindMolChiralCenters
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from multiprocessing import Pool, cpu_count
from itertools import combinations


def fused_ring_count(m):
    # count rings considering fused and spiro cycles as a single ring system
    # print(rings('C1CC23CCC2CC13'))  # 1
    # print(rings('O=C(N1CCNCC1)c1ccc(=O)oc1'))  # 2
    # print(rings('O=C(C1CCC(=O)C23CCCC2CCC13)N1CCNC2CCCC12'))  # 2
    # print(rings('O=C(C1CCC(=O)C23CCCC2CCC13)N1CCNC2C1CCC21CCCC1'))  # 2
    # print(rings('C1CC2(C1)CC1(C2)C2CCC22CCC12'))  # 1
    # print(rings('CC12CCC(C1)C1CCC21'))  # 1
    # print(rings('CC12CCC3(CCC3C1)C2'))  # 1
    # print(rings('CC'))  # 0
    # print(rings('C1CC2CCCC(C1)CCCC2'))  # 1
    q = m.GetRingInfo()
    rings = [set(r) for r in q.AtomRings()]
    go_next = True
    while go_next:
        go_next = False
        for i, j in combinations(range(len(rings)), 2):
            if rings[i] & rings[j]:
                q = rings[i] | rings[j]
                del rings[j], rings[i]
                rings.append(q)
                go_next = True
                break
    return len(rings)


def count_hbd_hba_atoms(m):
    HDonor = m.GetSubstructMatches(HDonorSmarts)
    HAcceptor = m.GetSubstructMatches(HAcceptorSmarts)
    return len(set(HDonor + HAcceptor))


def calc(smi, name):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        try:
            hba = rdMolDescriptors.CalcNumHBA(m)

            hbd = rdMolDescriptors.CalcNumHBD(m)
            nrings = rdMolDescriptors.CalcNumRings(m)
            rtb = rdMolDescriptors.CalcNumRotatableBonds(m)
            psa = rdMolDescriptors.CalcTPSA(m)
            logp, mr = rdMolDescriptors.CalcCrippenDescriptors(m)
            mw = rdMolDescriptors._CalcMolWt(m)
            csp3 = rdMolDescriptors.CalcFractionCSP3(m)
            hac = m.GetNumHeavyAtoms()
            if hac == 0:
                fmf = 0
            else:
                fmf = GetScaffoldForMol(m).GetNumHeavyAtoms() / hac
            qed = QED.qed(m)
            nrings_fused = fused_ring_count(m)
            n_unique_hba_hbd_atoms = count_hbd_hba_atoms(m)
            max_ring_size = len(max(m.GetRingInfo().AtomRings(), key=len, default=()))
            n_chiral_centers = len(FindMolChiralCenters(m, includeUnassigned=True))
            return name, hba, hbd, hba + hbd, nrings, rtb, round(psa, 2), round(logp, 2), round(mr, 2), round(mw, 2), \
                   round(csp3, 3), round(fmf, 3), round(qed, 3), hac, nrings_fused, n_unique_hba_hbd_atoms, \
                   max_ring_size, n_chiral_centers
        except:
            sys.stderr.write(f'molecule {name} was omitted due to an error in calculation of some descriptors\n')
            return None
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

    parser = argparse.ArgumentParser(description='Calculate some physicochemical parameters with RDKit.\n'
                                                 'HBA: number of H-bond acceptors\n'
                                                 'HBD: number of H-bond donors\n'
                                                 'complexity: HBA + HBD\n'
                                                 'NumRings: number of rings\n'
                                                 'RTB: number of rotatable bonds\n'
                                                 'TPSA: topological polar surface area\n'
                                                 'logP: lipophilicity\n'
                                                 'MR: molecular refraction\n'
                                                 'MW: molecular weight\n'
                                                 'Csp3: fraction of sp3 carbons\n'
                                                 'fmf: fraction of atoms belonging to Murcko framework\n'
                                                 'QED: quantitative estimate of drug-likeness\n'
                                                 'HAC: heavy atom count\n'
                                                 'NumRingsFused: number of rings considering fused and spirocycles as a single ring\n'
                                                 'unique_HBAD: number of unique H-bond acceptors and H-bond donors atoms\n'
                                                 'max_ring_size: maximum ring size in a molecule\n'
                                                 'ChiralCenters: number of chiral centers (assigned and unassigned)\n',
                                     formatter_class=RawTextHelpFormatter)
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

    HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
    HAcceptorSmarts = Chem.MolFromSmarts('[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),' +
                                         '$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),' +
                                         '$([nH0,o,s;+0])]')

    p = Pool(min(ncpu, cpu_count()))

    with open(out_fname, 'wt') as f:
        f.write('\t'.join(['Name', 'HBA', 'HBD', 'complexity', 'NumRings', 'RTB', 'TPSA', 'logP', 'MR', 'MW', 'Csp3',
                           'fmf', 'QED', 'HAC', 'NumRingsFused', 'unique_HBAD', 'max_ring_size',
                           'ChiralCenters']) + '\n')
        for i, res in enumerate(p.imap(calc_mp, read_smi(in_fname), chunksize=100)):
            if res:
                f.write('\t'.join(map(str, res)) + '\n')
            if verbose and i % 100 == 0:
                sys.stderr.write('\r%i molecules passed' % (i + 1))
                sys.stderr.flush()