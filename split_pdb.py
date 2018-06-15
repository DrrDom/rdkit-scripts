from rdkit import Chem
import argparse
import os


def main(pdb_fname):
    pdb = Chem.MolFromPDBFile(pdb_fname, sanitize=False)
    chains = Chem.SplitMolByPDBChainId(pdb)
    for name, mol in chains.items():
        with Chem.PDBWriter(os.path.splitext(pdb_fname)[0] + "_chain_%s.pdb" % name) as f:
            f.write(mol)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split PDB by chains and save to separate PDB files.')
    parser.add_argument('-i', '--in', metavar='input.pdb', required=True,
                        help='input PDB file.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": pdb_fname = v

    main(pdb_fname=pdb_fname)


