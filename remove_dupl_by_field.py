#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 10-06-2020
# copyright       : Pavel Polishchuk 2020
# license         : GPL3
#==============================================================================

import re
import argparse

patt = re.compile('> {1,2}<(.*)>( +\([0-9]+\))?')


def get_fields(molstr):
    # return dict of field names and values
    d = dict()
    i = 0
    while i < len(molstr):
        if patt.fullmatch(molstr[i].strip()):
            d[patt.sub('\\1', molstr[i].strip())] = molstr[i + 1].strip()
            i += 1
        i += 1
    return d


def main_params(in_fname, out_fname, title, field_name):

    values = set()

    with open(in_fname) as ifs, open(out_fname, 'wt') as ofs:

        molstr = []
        for line in ifs:
            molstr.append(line)
            if line.strip() == '$$$$':
                if title:
                    if molstr[0].strip() not in values:
                        values.add(molstr[0].strip())
                        ofs.writelines(molstr)
                else:
                    fields = get_fields(molstr)
                    if field_name in fields:
                        if fields[field_name] not in values:
                            values.add(fields[field_name])
                            ofs.writelines(molstr)
                    else:
                        print(f'missing field {field_name}')
                molstr = []


def main():
    parser = argparse.ArgumentParser(description='Remove entries with duplicated mol title or field value.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True,
                        help='output text file with values of specified fields')
    parser.add_argument('-t', '--title', action='store_true', default=False,
                        help='if set then molecule titles will be considered for searching of duplicates. '
                             'It has a priority over field_name. Only one line value will be considered.')
    parser.add_argument('-f', '--field_name', metavar='STRING', required=False, default=None,
                        help='field name to be used for searching of duplicates.')

    args = parser.parse_args()

    main_params(args.input, args.output, args.title, args.field_name)


if __name__ == '__main__':
    main()
