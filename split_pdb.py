from rdkit import Chem
import argparse
import os


def calc(pdb_fname):
    pdb = Chem.MolFromPDBFile(pdb_fname, sanitize=False)
    chains = Chem.SplitMolByPDBChainId(pdb)
    for name, mol in chains.items():
        w = Chem.PDBWriter(os.path.splitext(os.path.abspath(pdb_fname))[0] + "_chain_%s.pdb" % name)
        w.write(mol)
        w.close()


def main():
    parser = argparse.ArgumentParser(description='Split PDB by chains and save to separate PDB files.')
    parser.add_argument('-i', '--in', metavar='input.pdb', required=True,
                        help='input PDB file.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": pdb_fname = v

    calc(pdb_fname=pdb_fname)


if __name__ == '__main__':
    main()
