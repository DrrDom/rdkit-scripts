#!/usr/bin/env python
#==============================================================================
# author          : Pavel Polishchuk
# date            : 01-11-2014
# version         : 0.1
# python_version  : 3.2
# copyright       : Pavel Polishchuk 2014
# license         : GPL3
#==============================================================================

import re
import argparse
from itertools import chain

patt = re.compile('> {1,2}<(.*)>( +\([0-9]+\))?')


def get_fields(molstr):
    # return dict of field names and values
    d = dict()
    i = 0
    while i < len(molstr):
        if patt.fullmatch(molstr[i]):
            d[patt.sub('\\1', molstr[i])] = molstr[i + 1]
            i += 1
        i += 1
    return d


def main_params(in_fname, out_fname, title, field_names, all_fields, skip_value):

    with open(in_fname) as ifs:

        output = []

        molstr = []
        for line in ifs:
            l = line.strip()
            if l != '$$$$':
                molstr.append(l)
            else:
                d = dict()
                fields = get_fields(molstr)
                if title:
                    if molstr[0]:
                        d['Title'] = molstr[0]
                    else:
                        d['Title'] = ''
                if all_fields:
                    d.update(fields)
                else:
                    if field_names:
                        d.update({f: fields[f] for f in field_names if f in fields})
                output.append(d)
                molstr = []

    with open(out_fname, "w") as ofs:
        # get sorted unique field names
        field_names = sorted(list(set(list(chain.from_iterable([list(s.keys()) for s in output])))))
        if 'Title' in field_names:
            field_names.remove('Title')
            field_names.insert(0, 'Title')
        ofs.write("\t".join(field_names) + "\n")
        for item in output:
            line = [item.get(f, '') for f in field_names]
            if skip_value:
                # replace skipped values (NA) with empty string
                line = [v if v != skip_value else '' for v in line]
                # skip lines having all empty values
                if title:
                    if all(v == '' for v in line[1:]):
                        continue
                else:
                    if all(v == '' for v in line):
                        continue
            ofs.write("\t".join(line) + "\n")


def main():
    parser = argparse.ArgumentParser(description='Extract field values from sdf-files.')
    parser.add_argument('-i', '--in', metavar='FILENAME', required=True,
                        help='input SDF file, molecules should have titles')
    parser.add_argument('-o', '--out', metavar='FILENAME', required=True,
                        help='output text file with values of specified fields')
    parser.add_argument('-t', '--title', action='store_true', default=False,
                        help='If true then molecules titles will be extracted')
    parser.add_argument('-f', '--field_names', metavar='[field_name_1 field_name_2 ...]',
                        required=False, default=None, nargs='*',
                        help='space separated list of field names for extraction')
    parser.add_argument('-a', '--all_fields', action='store_true', default=False,
                        help='extract all fields.')
    parser.add_argument('--skip', default=None,
                        help='specify field value which will be skipped and replaced with empty string. '
                             'Usually NA can be set to skip missing values. '
                             'By default no values are skipped.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "title": title = v
        if o == "field_names": field_names = v
        if o == "all_fields": all_fields = v
        if o == "skip": skip_value = v

    main_params(in_fname, out_fname, title, field_names, all_fields, skip_value)


if __name__ == '__main__':
    main()
