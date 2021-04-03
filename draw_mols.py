#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import sys
from read_input import read_input
from rdkit import Chem
from multiprocessing import Pool, cpu_count

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import DrawingOptions
import pyvips

DrawingOptions.atomLabelFontSize = 300
DrawingOptions.elemDict = {}


def main(input_fname, output_fname, mols_per_file, mols_per_row, template):

    mols = []
    with open(input_fname) as f:
        for line in f:
            items = line.strip().split()
            m = Chem.MolFromSmiles(items[0])
            if m:
                mols.append([m] + items[1:])

    if template is not None:

        if template.endswith('.smi'):
            with open(template) as f:
                smi = f.readline().strip().split()[0]
                tmol = Chem.MolFromSmiles(ami)
                AllChem.Compute2DCoords(tmol)
        elif template.endswith('.mol'):
            tmol = Chem.MolFromMolFile(template)
        else:
            raise ValueError('reference file format is wrong. Only SMI or MOL files are allowed.')

        # check matching mols if reference was given
        matched_mols = []
        for items in mols:
            if not items[0].HasSubstructMatch(tmol):
                print(f'{items[1]} does not match template and will be omitted')
            else:
                matched_mols.append(items)
        mols = matched_mols

        for items in mols:
            AllChem.GenerateDepictionMatching2DStructure(items[0], tmol)

    if output_fname.endswith('.png'):
        output_fname = output_fname[:-4]

    if mols_per_file is None:
        mols_per_file = len(mols)

    for j, i in enumerate(range(0, len(mols), mols_per_file)):
        selected_mols = [items[0] for items in mols[i:i + mols_per_file]]
        legend = []
        for items in mols[i:i + mols_per_file]:
            legend.append(' / '.join(items[1:3]))
        img = Draw.MolsToGridImage(selected_mols, molsPerRow=mols_per_row, subImgSize=(300, 200), legends=legend,
                                   useSVG=True)
        # with open(output_fname + str(j).zfill(3) + '.png', 'w') as f:
        #     f.write(img)
        pyvips.Image.svgload_buffer(img.encode('utf-8'), dpi=300).write_to_file(output_fname + str(j).zfill(3) + '.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create images from input molecules. If a reference structure was '
                                                 'supplied all molecule coordinates will be aligned correspondingly.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='input SMILES file. No header. Separator whitespaces.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=str,
                        help='output PNG file.')
    parser.add_argument('-m', '--mols_per_file', metavar='INTEGER', required=False, default=None, type=int,
                        help='number of molecules per file.')
    parser.add_argument('-r', '--mols_per_row', metavar='INTEGER', required=False, default=5, type=int,
                        help='number of molecules per row. Default: 5.')
    parser.add_argument('-t', '--template', metavar='FILENAME', required=False, default=None, type=str,
                        help='MOL or SMI file containing a reference molecule used as template for coordinate '
                             'alignment.')
    args = parser.parse_args()

    main(args.input, args.output, args.mols_per_file, args.mols_per_row, args.template)
