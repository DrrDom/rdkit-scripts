#!/usr/bin/env python

import argparse
import os
import sqlite3


def main():
    parser = argparse.ArgumentParser(description='Extract mol blocks of specified mol ids and additional fields into '
                                                 'SDF. Also it is possible to extract SMILES file.')
    parser.add_argument('-i', '--input', metavar='input.db', required=True, type=str,
                        help='SQLite DB, which is output of vina_dock script.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True, type=str,
                        help='output SDF file (with mol blocks) or SMILES. Output format is guessed from extension.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None,
                        help='comma separated list of mol ids in DB or a text file with mol ids on individual lines. '
                             'If omitted all records in DB will be saved to SDF.')
    parser.add_argument('-f', '--first_entry', action='store_true', default=False,
                        help='retrieve only the first entry of each molecule from the database.')
    parser.add_argument('--add_sql', default=None, type=str,
                        help='sql string which will be added to the sql query to select data. This may be useful '
                             'to make additonal filtering of database compounds, e.g. iteration > 0.')
    parser.add_argument('--fields', default=[], type=str, nargs="*",
                        help='names of fields in database to additionally retrieve.')

    args = parser.parse_args()

    if args.ids is None:
        ids = tuple()
    elif os.path.isfile(args.ids):
        with open(args.ids) as f:
            ids = [line.strip() for line in f]
    else:
        ids = args.ids.split(',')

    conn = sqlite3.connect(args.input)
    cur = conn.cursor()

    # True - if the table tautomers exists and it has at least one row
    try:
        tautomers_exist = cur.execute('SELECT COUNT(rowid) FROM tautomers').fetchone()[0] > 0
    except sqlite3.Error:
        tautomers_exist = False

    ext = args.output.rsplit('.', 1)[1].lower()
    if ext == 'sdf':
        main_field = 'mol_block'
    elif ext == 'smi':
        main_field = 'smi'
    else:
        raise ValueError('Wrong extension of output file. Only SDF and SMI are allowed.')

    if args.fields:
        sql = f"SELECT {main_field}, {','.join(args.fields)} FROM mols WHERE mol_block IS NOT NULL"
    else:
        sql = f"SELECT {main_field} FROM mols WHERE mol_block IS NOT NULL"
    if args.ids:
        sql += f" AND id IN ({','.join('?' * len(ids))})"
    if args.add_sql:  # for example: iteration > 0
        sql += f" AND {args.add_sql}"
    if tautomers_exist:
        sql += f" AND id NOT IN (SELECT id FROM tautomers)"

    if tautomers_exist:
        sql += " UNION "
        if args.fields:
            sql += f" SELECT {main_field}, {','.join(args.fields)} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        else:
            sql += f" SELECT {main_field} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        if args.ids:
            sql += f" AND id IN ({','.join('?' * len(ids))})"

    if args.first_entry:
        sql += " GROUP BY id HAVING MIN(rowid)"

    if tautomers_exist:
        res = cur.execute(sql, ids + ids)   # ids should be duplicated to be selected from mols and tautomers tables
    else:
        res = cur.execute(sql, ids)  # ids can be empty if we retrieve all molecules, this works

    with open(args.output, 'wt')as f:
        for item in res:  # (id, mol_block, ...)
            if ext == 'sdf':
                mol_block = item[0]
                f.write(mol_block)
                for prop_name, prop_value in zip(args.fields, item[1:]):
                    f.write(f'>  <{prop_name}>\n')
                    f.write(f'{str(prop_value)}\n\n')
                f.write('$$$$\n')
            elif ext == 'smi':
                f.write('\t'.join(map(str, item)) + '\n')


if __name__ == '__main__':
    main()