#!/usr/bin/env python

import argparse
import os
import sqlite3


def main():
    parser = argparse.ArgumentParser(description='Extract mol blocks of specified mol ids into SDF file.')
    parser.add_argument('-i', '--input', metavar='input.db', required=True, type=str,
                        help='SQLite DB, which is output of vina_dock script.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True, type=str,
                        help='output SDF file.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None,
                        help='comma separated list of mol ids in DB or a text file with mol ids on individual lines. '
                             'If omitted all records in DB will be saved to SDF.')

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
    sql = "SELECT mol_block, docking_score FROM mols"
    if ids is not None:
        sql += f" WHERE id IN ({','.join(['?'] * len(ids))})"
        res = cur.execute(sql, ids)
    else:
        res = cur.execute(sql)

    with open(args.output, 'wt')as f:
        for mol_block, docking_score in res:
            f.write(mol_block)
            f.write('>  <docking_score>\n')
            f.write(f'{docking_score}\n\n')
            f.write('$$$$\n')


if __name__ == '__main__':
    main()
