import argparse
import os
import sqlite3
import sys
import tempfile
from functools import partial
from multiprocessing import Pool, Manager

import dask
import shutil
import random, string
from os import system
from dask import bag
from dask.distributed import Lock as daskLock, Client
from preparation_for_docking import create_db, save_sdf, add_protonation, ligand_preparation, \
    fix_pdbqt, pdbqt2molblock, cpu_type, filepath_type


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def save_pdbqt_to_file(pdbqt, mol_id, dir_name):
    ligand_pdbqt_file = os.path.join(dir_name, f'{mol_id}.pdbqt')
    with open(ligand_pdbqt_file, 'wt') as f:
        f.write(pdbqt)
    return ligand_pdbqt_file


def docking(script_file, ligand_pdbqt_file, ligand_out_fname, receptor_pdbqt_fname, protein_setup, exhaustiveness, seed,
            scoring, addH, cnn_scoring, cnn, ncpu):

    system(
        f'{script_file} --receptor {receptor_pdbqt_fname} --ligand {ligand_pdbqt_file} --out {ligand_out_fname} '
        f'--config {protein_setup} --exhaustiveness {exhaustiveness} --seed {seed} --scoring {scoring} --cpu {ncpu} '
        f'--addH {addH} --cnn_scoring {cnn_scoring} --cnn {cnn}')


def get_pdbqt_and_score(ligand_out_fname):
    pdbqt_out = open(ligand_out_fname).read()
    score = round(float(pdbqt_out.split(' ')[3].split('REMARK')[0]), 3)
    return score, pdbqt_out


def process_mol_docking(mol_id, smi, script_file, tmpdir, receptor_pdbqt_fname, protein_setup, dbname, seed, exhaustiveness,
                        scoring, addH, cnn_scoring, cnn, ncpu, lock=None):

    def insert_data(dbname, pdbqt_out, score, mol_block, mol_id):
        with sqlite3.connect(dbname) as conn:
            conn.execute("""UPDATE mols
                               SET pdb_block = ?,
                                   docking_score = ?,
                                   mol_block = ?,
                                   time = CURRENT_TIMESTAMP
                               WHERE
                                   id = ?
                            """, (pdbqt_out, score, mol_block, mol_id))

    ligand_pdbqt = ligand_preparation(smi, seed)
    if ligand_pdbqt is None:
        return mol_id

    ligand_pdbqt_file = save_pdbqt_to_file(ligand_pdbqt, mol_id, tmpdir)
    ligand_out_fname = ligand_pdbqt_file.rsplit('.', 1)[0] + '_dock.pdbqt'  # saving result after docking to file
    print('ligand_out_fname', ligand_out_fname)

    docking(script_file=script_file, ligand_pdbqt_file=ligand_pdbqt_file, ligand_out_fname=ligand_out_fname,
                               receptor_pdbqt_fname=receptor_pdbqt_fname, protein_setup=protein_setup,
                               exhaustiveness=exhaustiveness, seed=seed, scoring=scoring, addH=addH,
                               cnn_scoring=cnn_scoring, cnn=cnn, ncpu=ncpu)

    score, pdbqt_out = get_pdbqt_and_score(ligand_out_fname)

    mol_block = pdbqt2molblock(pdbqt_out, smi, mol_id)
    if mol_block is None:
        pdbqt_out = fix_pdbqt(pdbqt_out)
        mol_block = pdbqt2molblock(pdbqt_out, smi, mol_id)
        if mol_block:
            sys.stderr.write('PDBQT was fixed\n')

    if lock is not None:  # multiprocessing
        with lock:
            insert_data(dbname, pdbqt_out, score, mol_block, mol_id)
    else:  # dask
        with daskLock(dbname):
            insert_data(dbname, pdbqt_out, score, mol_block, mol_id)

    return mol_id


def iter_docking(script_file, tmpdir, dbname, receptor_pdbqt_fname, protein_setup, protonation, exhaustiveness, seed,
                 scoring, addH, cnn_scoring, cnn, ncpu, use_dask, add_sql=None):
    '''
    This function should update output db with docked poses and scores. Docked poses should be stored as pdbqt (source)
    and mol block. All other post-processing will be performed separately.
    :param script_file: path to gnina file
    :param tmpdir: dir, where temporary files (such as .pdbqt and _dock.pdbqt) will be stored
    :param dbname:
    :param receptor_pdbqt_fname:
    :param protein_setup: text file with vina grid box parameters
    :param protonation: True or False
    :param exhaustiveness: int
    :param seed: int
    :param scoring: type of scoring, for example: 'vina'
    :param ncpu: int
    :param use_dask: indicate whether or not using dask cluster
    :type use_dask: bool
    :param add_sql: string with additional selection requirements which will be concatenated to the main SQL query
                    with AND operator, e.g. "iteration = 1".
    :return:
    '''

    with sqlite3.connect(dbname) as conn:
        cur = conn.cursor()
        smi_field_name = 'smi_protonated' if protonation else 'smi'
        sql = f"SELECT id, {smi_field_name} FROM mols WHERE docking_score IS NULL"
        if isinstance(add_sql, str) and add_sql:
            sql += f" AND {add_sql}"
        smiles_dict = dict(cur.execute(sql))
    if not smiles_dict:
        return

    if use_dask:
        i = 0
        b = bag.from_sequence(smiles_dict.items(), npartitions=len(smiles_dict))
        for i, mol_id in enumerate(b.starmap(process_mol_docking, script_file=script_file, tmpdir=tmpdir,
                                             receptor_pdbqt_fname=receptor_pdbqt_fname, protein_setup=protein_setup,
                                             dbname=dbname, exhaustiveness=exhaustiveness,
                                             seed=seed, scoring=scoring, addH=addH, ncpu=1).compute(),
                                   1):
            if i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were docked')
        sys.stderr.write(f'\r{i} molecules were docked\n')

    else:
        pool = Pool(ncpu)
        manager = Manager()
        lock = manager.Lock()
        i = 0
        for i, mol_id in enumerate(pool.starmap(partial(process_mol_docking, script_file=script_file, tmpdir=tmpdir, dbname=dbname,
                                                        receptor_pdbqt_fname=receptor_pdbqt_fname, protein_setup=protein_setup,
                                                        seed=seed, exhaustiveness=exhaustiveness, scoring=scoring, addH=addH,
                                                        cnn_scoring=cnn_scoring, cnn=cnn, ncpu=1, lock=lock),
                                                smiles_dict.items()), 1):
            if i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were docked')
        sys.stderr.write(f'\r{i} molecules were docked\n')


def main():
    parser = argparse.ArgumentParser(description='Perform docking of input molecules using GNINA. The script '
                                                 'automates the whole pipeline: protonate molecules, creates '
                                                 '3D structures, converts to PDBQT format, run docking using a single '
                                                 'machine (multiprocessing) or a cluster of servers (dask), stores '
                                                 'the best scores and poses in PDBQT and MOL formats to DB.\n\n'
                                                 'It has multiple dependencies:\n'
                                                 '  - rdkit - conda install -c conda-forge rdkit\n'
                                                 '  - gnina - wget https://github.com/gnina/gnina/releases/download/v1.0.1/gnina\n'
                                                 'possible problems during installation, for example, version "CXXABI_1.3.8" not found (required by ./gnina)\n'
                                                 'you need install libstdc++.so.6 and add way to this file in .bashrc\n'
                                                 'for example: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/lib64/\n'
                                                 '  - openbabel - conda install -c conda-forge openbabel\n'
                                                 '  - meeko - pip install git+https://github.com/forlilab/Meeko@7b1a60d9451eabaeb16b08a4a497cf8e695acc63\n'
                                                 '  - Chemaxon cxcalc utility\n\n'
                                                 'To run on a single machine:\n'
                                                 '  gnina_docking.py -i input.smi -o output.db -p protein.pdbqt -s config.txt -c 4 -v\n\n'
                                                 'To run on several machines using dask ssh-cluster (on PBS system):\n'
                                                 '  dask-ssh --hostfile $PBS_NODEFILE --nprocs 32 --nthreads 1 &\n'
                                                 '  sleep 10\n'
                                                 '  gnina_docking.py -i input.smi -o output.db -p protein.pdbqt -s config.txt -c 4 -v --hostfile $PBS_NODEFILE\n\n'
                                                 '  config.txt contains coordinates of a gridbox\n'
                                                 '  $PBS_NODEFILE contains the list of addresses of servers\n',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with molecules (SMI, SDF, SDF.GZ, PKL). Maybe be omitted if output DB exists. '
                             'In this case calculations will be continued from interrupted point.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with scores and poses in PDBQT and MOL formats. It also stores '
                             'other information (input structures, protein pdbqt file and grid box config). '
                             'If output DB exists all other inputs will be ignored and calculations will be continued.')
    parser.add_argument('--no_protonation', action='store_true', default=False,
                        help='disable protonation of molecules before docking. Protonation requires installed '
                             'cxcalc chemaxon utility. It will be omitted if output DB exists.')
    parser.add_argument('-p', '--protein', metavar='protein.pdbqt', required=False, type=filepath_type,
                        help='input PDBQT file with a prepared protein. It will be omitted if output DB exists.')
    parser.add_argument('-s', '--protein_setup', metavar='protein.log', required=False, type=filepath_type,
                        help='input text file with GNINA docking setup. It will be omitted if output DB exists.')
    parser.add_argument('-e', '--exhaustiveness', metavar='INTEGER', required=False, type=int, default=8,
                        help='exhaustiveness of docking search.')
    parser.add_argument('--seed', metavar='INTEGER', required=False, type=int, default=0,
                        help='seed to make results reproducible.')
    parser.add_argument('--sdf', action='store_true', default=False,
                        help='save best docked poses to SDF file with the same name as output DB. Can be used with DB '
                             'of previously docked molecules to retrieve SDF file.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=filepath_type, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('--tmpdir', metavar='DIRNAME', required=False, type=filepath_type, default=None,
                        help='path to a dir where to store temporary files accessible to a program. If use dask this '
                             'argument must be specified because dask cannot access ordinary tmp locations.')
    parser.add_argument('--prefix', metavar='STRING', required=False, type=str, default=None,
                        help='prefix which will be added to all names. This might be useful if multiple runs are made '
                             'which will be analyzed together.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    parser.add_argument('--install_dir', metavar='DIRNAME', type=filepath_type, required=True,
                        help='path to the dir with installed GNINA.')
    parser.add_argument('--scoring', metavar='STRING', required=True,
                        choices=['ad4_scoring', 'default', 'dkoes_fast', 'dkoes_scoring', 'dkoes_scoring_old', 'vina', 'vinardo'],
                        help='type of scoring function.')
    parser.add_argument('--addH', action='store_true', default=False,
                        help='add hydrogens in ligands.')
    parser.add_argument('--cnn_scoring', metavar='STRING', required=False, default='rescore', choices=['None', 'rescore', 'refinement', 'all'],
                        help='type of CNN scoring.')
    parser.add_argument('--cnn', metavar='STRING', required=False, default='dense_ensemble', choices=['crossdock_default2018_ensemble', 'dense_ensemble'],
                        help='type of built-in model to use.')


    args = parser.parse_args()


    gnina_script_dir = os.path.join(args.install_dir, 'gnina')

    if args.tmpdir is None:
        tmpdir = (os.path.join(os.path.dirname(args.output),''.join(random.sample(string.ascii_lowercase, 6))))
        os.makedirs(tmpdir, exist_ok=True)
    else:
        tmpdir = args.tmpdir
        tempfile.tempdir = args.tmpdir

    if args.hostfile is not None:
        dask.config.set({'distributed.scheduler.allowed-failures': 30})
        dask_client = Client(open(args.hostfile).readline().strip() + ':8786')

    if not os.path.isfile(args.output):
        create_db(args.output, args.input, not args.no_protonation, args.protein, args.protein_setup, args.prefix)

    add_protonation(args.output)

    conn = sqlite3.connect(args.output)
    protein = tempfile.NamedTemporaryFile(suffix='.pdbqt', mode='w', encoding='utf-8')
    protein.write(list(conn.execute('SELECT protein_pdbqt FROM setup'))[0][0])
    protein.flush()
    setup = tempfile.NamedTemporaryFile(suffix='.txt', mode='w', encoding='utf-8')
    setup.write(list(conn.execute('SELECT protein_setup FROM setup'))[0][0])
    setup.flush()
    protonation = list(conn.execute('SELECT protonation FROM setup'))[0][0]

    iter_docking(script_file=gnina_script_dir, tmpdir=tmpdir, dbname=args.output, receptor_pdbqt_fname=protein.name, protein_setup=setup.name,
                 protonation=protonation, exhaustiveness=args.exhaustiveness, seed=args.seed, scoring=args.scoring, addH=args.addH,
                 cnn_scoring=args.cnn_scoring, cnn=args.cnn, ncpu=args.ncpu, use_dask=args.hostfile is not None)

    if args.sdf:
        save_sdf(args.output)

    if args.tmpdir is None:
        shutil.rmtree(tmpdir, ignore_errors=True) # to delete temporary dir with files adter docking


if __name__ == '__main__':
    main()