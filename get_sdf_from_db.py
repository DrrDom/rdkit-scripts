#!/usr/bin/env python

import argparse
import os
import sqlite3


def main():
    parser = argparse.ArgumentParser(description='Extract mol blocks of specified mol ids into SDF file and extract '
                                                 'their docking scores.')
    parser.add_argument('-i', '--input', metavar='input.db', required=True, type=str,
                        help='SQLite DB, which is output of vina_dock script.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True, type=str,
                        help='output SDF file.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None,
                        help='comma separated list of mol ids in DB or a text file with mol ids on individual lines. '
                             'If omitted all records in DB will be saved to SDF.')
    parser.add_argument('-f', '--first_entry', action='store_true', default=False,
                        help='retrieve only the first entry of each molecule from the database.')
    parser.add_argument('-x', '--no_fields', action='store_true', default=False,
                        help='choose this option if you do not want to retrieve any further fields from a database.')

    args = parser.parse_args()

    if args.ids is None:
        ids = None
    elif os.path.isfile(args.ids):
        with open(args.ids) as f:
            ids = [line.strip() for line in f]
    else:
        ids = args.ids.split(',')

    conn = sqlite3.connect(args.input)
    cur = conn.cursor()
    if args.no_fields:
        sql = "SELECT mol_block FROM mols"
    else:
        sql = "SELECT mol_block, docking_score FROM mols"
    if ids is not None:
        sql += f" WHERE id IN ({','.join(['?'] * len(ids))})"
    if args.first_entry:
        sql += " GROUP BY id HAVING MIN(rowid) ORDER BY rowid"

    if ids is not None:
        res = cur.execute(sql, ids)
    else:
        res = cur.execute(sql)

    with open(args.output, 'wt')as f:
        for item in res:
            mol_block = item[0]
            f.write(mol_block)
            if len(item) == 2:
                docking_score = item[1]
                f.write('>  <docking_score>\n')
                f.write(f'{docking_score}\n\n')
            f.write('$$$$\n')


if __name__ == '__main__':
    main()
