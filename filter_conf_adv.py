#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import numpy as np
import sys
from functools import partial
from read_input import read_input, assign_conf_props_to_mol
from rdkit import Chem
from rdkit.Chem.AllChem import GetConformerRMSMatrix
from multiprocessing import Pool, cpu_count
from sklearn.cluster import AgglomerativeClustering


def remove_confs_rms(mol, rms=0.25, keep_nconf=None, preferably_keep=None):
    """
    The function uses AgglomerativeClustering to select conformers.

    :param mol: input molecule with multiple conformers
    :param rms: discard conformers which are closer than given value to a kept conformer
    :param keep_nconf: keep at most the given number of conformers. This parameter has precedence over rms
    :param preferably_keep: a tuple with the field name and the value to preferably keep conformers having this value
    :return:
    """

    def gen_ids(ids):
        for i in range(1, len(ids)):
            for j in range(0, i):
                yield j, i

    if keep_nconf and mol.GetNumConformers() <= keep_nconf:
        return mol

    if mol.GetNumConformers() <= 1:
        return mol

    mol_tmp = Chem.RemoveHs(mol)   # calc rms for heavy atoms only
    rms_ = GetConformerRMSMatrix(mol_tmp, prealigned=True)

    cids = [c.GetId() for c in mol_tmp.GetConformers()]
    arr = np.zeros((len(cids), len(cids)))
    for (i, j), v in zip(gen_ids(cids), rms_):
        arr[i, j] = v
        arr[j, i] = v
    if keep_nconf:
        cl = AgglomerativeClustering(n_clusters=keep_nconf, linkage='complete', metric='precomputed').fit(arr)
    else:
        cl = AgglomerativeClustering(n_clusters=None, linkage='complete', metric='precomputed', distance_threshold=rms).fit(arr)

    keep_cids = []
    for i in set(cl.labels_):
        ids = np.where(cl.labels_ == i)[0]
        if preferably_keep is not None:
            ids_keep = []  # these are sequential ids in ids var
            for k, conf_id in enumerate(ids):
                try:
                    if mol_tmp.GetConformer(int(conf_id)).GetProp(preferably_keep[0]) == preferably_keep[1]:
                        ids_keep.append(k)
                except KeyError:
                    print(f'{mol_tmp.GetProp("_Name")}  {conf_id}')
                    print(f'\n{mol_tmp.GetConformer(int(conf_id)).GetPropsAsDict()}')
            if ids_keep:
                j = arr[np.ix_(ids, ids)].mean(axis=0)[ids_keep].argmin()
                j = ids_keep[j]
            else:
                j = arr[np.ix_(ids, ids)].mean(axis=0).argmin()
        else:
            j = arr[np.ix_(ids, ids)].mean(axis=0).argmin()
        keep_cids.append(cids[ids[j]])
    remove_ids = set(cids) - set(keep_cids)

    for cid in sorted(remove_ids, reverse=True):
        mol_tmp.RemoveConformer(cid)

    return mol_tmp


def calc(input_fname, output_fname, rms, nconf, preferably_keep, ncpu, verbose):

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    pool = Pool(max(min(cpu_count(), ncpu), 1))

    w = Chem.SDWriter(output_fname)
    for i, mol in enumerate(pool.imap(partial(remove_confs_rms,
                                              rms=rms,
                                              keep_nconf=nconf,
                                              preferably_keep=preferably_keep),
                                      (mol for mol, mol_id in read_input(input_fname, sdf_confs=True))),
                            1):
        for conf in mol.GetConformers():
            assign_conf_props_to_mol(conf, mol)
            w.write(mol, confId=conf.GetId())
        if verbose:
            sys.stderr.write(f'\rMolecules passed: {i}')
    if verbose:
        sys.stderr.write('\n')


def main():
    parser = argparse.ArgumentParser(description='Filter out conformers in the input file. Conformers are clustered '
                                                 'by RMSD and from each cluster a representative conformer is '
                                                 'selected. Selected conformers can be biased by specifying '
                                                 'a particular value of a conformer property (data field).')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SDF.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output SDF file.')
    parser.add_argument('-r', '--rms', metavar='NUMERIC', default=1, type=float,
                        help='RMS threshold to discard conformers. Default: 1 (Angstrom).')
    parser.add_argument('-n', '--nconf', metavar='INTEGER', default=None, type=int,
                        help='number of selected conformers. It has a precedence over an RMS argument. '
                             'Default: None (no filtering).')
    parser.add_argument('-k', '--preferably_keep', metavar=('STRING', 'STRING'), default=None, nargs=2,
                        help='name of the property field and its value to use to preferably keep corresponding '
                             'conformers having this property value.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=int,
                        help='number of cpu to use for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    args = parser.parse_args()

    calc(args.input, args.output, args.rms, args.nconf, args.preferably_keep, args.ncpu, args.verbose)


if __name__ == '__main__':
    main()
