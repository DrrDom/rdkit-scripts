import argparse
import os
import sys
import string
import random
import itertools
from multiprocessing import Pool, cpu_count
from functools import partial
# import numpy as np
from collections import defaultdict, Counter

from read_input import read_input
from pmapper.utils import load_multi_conf_mol


class SvmSaver:

    def __init__(self, file_name):
        self.__fname = file_name
        self.__varnames_fname = os.path.splitext(file_name)[0] + '.colnames'
        self.__molnames_fname = os.path.splitext(file_name)[0] + '.rownames'
        self.__varnames = dict()
        if os.path.isfile(self.__fname):
            os.remove(self.__fname)
        if os.path.isfile(self.__molnames_fname):
            os.remove(self.__molnames_fname)
        if os.path.isfile(self.__varnames_fname):
            os.remove(self.__varnames_fname)

    def save_mol_descriptors(self, mol_name, mol_descr_dict):

        new_varnames = list(mol_descr_dict.keys() - self.__varnames.keys())
        for v in new_varnames:
            self.__varnames[v] = len(self.__varnames)

        values = {self.__varnames[k]: v for k, v in mol_descr_dict.items()}

        if values:  # values can be empty if all descriptors are zero

            with open(self.__molnames_fname, 'at') as f:
                f.write(mol_name + '\n')

            if new_varnames:
                with open(self.__varnames_fname, 'at') as f:
                    f.write('\n'.join(new_varnames) + '\n')

            with open(self.__fname, 'at') as f:
                values = sorted(values.items())
                values_str = ('%i:%i' % (i, v) for i, v in values)
                f.write(' '.join(values_str) + '\n')

            return tuple(i for i, v in values)

        return tuple()


# class SvmSaver2:
#
#     def __init__(self, file_name):
#         self.__fname = file_name
#         self.__varnames_fname = os.path.splitext(file_name)[0] + '.colnames'
#         self.__molnames_fname = os.path.splitext(file_name)[0] + '.rownames'
#         self.__varnames = np.array([])
#         if os.path.isfile(self.__fname):
#             os.remove(self.__fname)
#         if os.path.isfile(self.__molnames_fname):
#             os.remove(self.__molnames_fname)
#         if os.path.isfile(self.__varnames_fname):
#             os.remove(self.__varnames_fname)
#
#     def save_mol_descriptors(self, mol_name, mol_descr_dict):
#
#         if mol_descr_dict:
#
#             names, values = zip(*sorted(mol_descr_dict.items()))
#             names = np.array(names)
#             a = names[~np.isin(names, self.__varnames)]
#             if a.size:
#                 self.__varnames = np.concatenate([self.__varnames, a])
#
#             sorter = np.argsort(self.__varnames)
#             ids = sorter[np.searchsorted(self.__varnames, names, sorter=sorter)]
#
#             with open(self.__molnames_fname, 'at') as f:
#                 f.write(mol_name + '\n')
#
#             if a.size:
#                 with open(self.__varnames_fname, 'at') as f:
#                     f.write('\n'.join(a) + '\n')
#
#             with open(self.__fname, 'at') as f:
#                 f.write(' '.join(f'{i}:{j}' for i, j in sorted(zip(ids, values))) + '\n')


# svm = SvmSaver2('test/1.txt')
#
# svm.save_mol_descriptors('mol1', {'d1': 3, 'd3': 5, 'd8': 1})
# svm.save_mol_descriptors('mol2', {'d2': 1, 'd3': 2})


def process_mol(mol, mol_title, descr_num):
    # descr_num - list of int
    ps = load_multi_conf_mol(mol)
    res = []
    for p in ps:
        tmp = dict()
        for n in descr_num:
            tmp.update(p.get_descriptors(ncomb=n))
        res.append(tmp)
    ids = [c.GetId() for c in mol.GetConformers()]
    ids, res = zip(*sorted(zip(ids, res)))  # reorder output by conf ids
    return mol_title, res


def process_mol_map(items, descr_num):
    return process_mol(*items, descr_num=descr_num)


def main():
    parser = argparse.ArgumentParser(description='Calculate 3D pharmacophore descriptors and remove rarely '
                                                 'occurred ones. Descriptors are generated using binning step 1A. '
                                                 'A temporary file is created containing all descriptors which are '
                                                 'filtered to create an output file. In the case of PKL input file '
                                                 'the output order of conformers will be according to conformer ids '
                                                 'in an RDKit Mol object (not their actual order).',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='SDF or PKL file. In the case of SDF molecular titles will be used to identify molecular '
                             'instances. PKL file should contain tuples of molecules and their titles. SDF file is not '
                             'implemented yet.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='text file with a computed descriptor matrix.')
    parser.add_argument('-d', '--descr', metavar='INTEGER', default=[4], nargs='+', type=int,
                        help='number of features in a single descriptor. Can be set from 1 to 4. Multiple entries are '
                             'allowed. Default: 4.')
    parser.add_argument('-r', '--remove', metavar='NUMERIC', required=False, default=0.05, type=float,
                        help='minimal percentage of compounds with non-zero descriptor values to keep this descriptor '
                             'in the output matrix. Default: 0.05 (means 5%).')
    parser.add_argument('-t', '--keep_temp', action='store_true', default=False,
                        help='whether to not remove temporary files with descriptors containing all descriptors.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1, type=int,
                        help='number of cores for calculation. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    if args.remove < 0 or args.remove > 1:
        raise ValueError('Value of the "remove" argument is out of range [0, 1]')

    for v in args.descr:
        if v < 1 or v > 4:
            raise ValueError('The number of features in a single descriptor should be within 1-4 range.')

    pool = Pool(max(min(args.ncpu, cpu_count()), 1))

    tmp_fname = os.path.splitext(args.output)[0] + '.' + ''.join(random.sample(string.ascii_lowercase, 6)) + '.txt'
    svm = SvmSaver(tmp_fname)

    stat = defaultdict(set)

    # create temp file with all descriptors
    for i, (mol_title, desc) in enumerate(pool.imap(partial(process_mol_map, descr_num=args.descr),
                                                    read_input(args.input), chunksize=1), 1):
        print(mol_title, len(desc))
        for desc_dict in desc:
            if desc_dict:
                ids = svm.save_mol_descriptors(mol_title, desc_dict)
                stat[mol_title].update(ids)
        if args.verbose and i % 10 == 0:
            sys.stderr.write(f'\r{i} molecule records were processed')
    sys.stderr.write('\n')

    if args.remove == 0:  # if no remove - rename temp files to output files
        os.rename(tmp_fname, args.output)
        os.rename(os.path.splitext(tmp_fname)[0] + '.colnames', os.path.splitext(args.output)[0] + '.colnames')
        os.rename(os.path.splitext(tmp_fname)[0] + '.rownames', os.path.splitext(args.output)[0] + '.rownames')

    else:
        # determine frequency of descriptors occurrence and select frequently occurred
        c = Counter(itertools.chain.from_iterable(stat.values()))
        threshold = len(stat) * args.remove
        desc_ids = {k for k, v in c.items() if v >= threshold}

        # create output files with removed descriptors

        replace_dict = dict()  # old_id, new_id
        with open(os.path.splitext(args.output)[0] + '.colnames', 'wt') as fout:
            with open(os.path.splitext(tmp_fname)[0] + '.colnames') as fin:
                for i, line in enumerate(fin):
                    if i in desc_ids:
                        replace_dict[i] = len(replace_dict)
                        fout.write(line)

        with open(os.path.splitext(args.output)[0] + '.rownames', 'wt') as fmol, open(args.output, 'wt') as ftxt:
            with open(os.path.splitext(tmp_fname)[0] + '.rownames') as fmol_tmp, open(tmp_fname) as ftxt_tmp:
                for line1, line2 in zip(fmol_tmp, ftxt_tmp):
                    desc_str = []
                    for item in line2.strip().split(' '):
                        i, v = item.split(':')
                        i = int(i)
                        if i in replace_dict:
                            desc_str.append(f'{replace_dict[i]}:{v}')
                    if desc_str:
                        fmol.write(line1)
                        ftxt.write(' '.join(desc_str) + '\n')

        if not args.keep_temp:
            os.remove(tmp_fname)
            os.remove(os.path.splitext(tmp_fname)[0] + '.colnames')
            os.remove(os.path.splitext(tmp_fname)[0] + '.rownames')


if __name__ == '__main__':
    main()
