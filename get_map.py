#!/usr/bin/env python3
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from multiprocessing import cpu_count
import numpy as np
import umap
from openTSNE import TSNE
from read_input import read_input


def compute_fps(mols, fp_type="morgan", nBits=2048, radius=2):
    fps = []
    for mol in mols:
        if mol is None:
            fps.append(np.zeros((nBits,), dtype=int))
            continue
        if fp_type == "morgan":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        elif fp_type == "atom_pair":
            fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nBits)
        elif fp_type == "torsion":
            fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nBits)
        elif fp_type == "rdkit":
            fp = Chem.RDKFingerprint(mol, fpSize=nBits)
        else:
            raise ValueError(f"Unsupported fingerprint type: {fp_type}")
        arr = np.zeros((1,), dtype=int)
        Chem.DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr)
    return np.array(fps)


def main():
    parser = argparse.ArgumentParser(description="2D UMAP projection of molecules")
    parser.add_argument("-i", "--input", required=True, help="Input SMILES file")
    parser.add_argument("-o", "--output", required=True, help="Output coordinates file")
    parser.add_argument("-m", "--method", required=False, choices=['umap', 'tsne'], default='umap',
                        help="Output coordinates file")
    parser.add_argument("-f", "--fingerprint", choices=["morgan", "atom_pair", "torsion", "rdkit"],
                        default="morgan", help="Type of fingerprint to compute")
    parser.add_argument("--radius", type=int, default=2, help="Radius for Morgan fingerprints")
    parser.add_argument("--n_bits", type=int, default=2048, help="Number of bits in fingerprint")
    parser.add_argument("--n_neighbors", type=float, default=15, help="UMAP n_neighbors / t-SNE perplexity")
    parser.add_argument("--min_dist", type=float, default=0.1, help="UMAP min_dist parameter")
    parser.add_argument("--metric", type=str, default="jaccard", help="UMAP distance metric")
    args = parser.parse_args()

    mols, mol_names = zip(*(read_input(args.input)))
    fps = compute_fps(mols, args.fingerprint, args.n_bits, args.radius)

    if args.method == "umap":
        reducer = umap.UMAP(
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            metric=args.metric,
            n_components=2
        )
        coords = reducer.fit_transform(fps)

    elif args.method == "tsne":
        reducer = TSNE(n_jobs=cpu_count(),
                       perplexity=args.n_neighbors,
                       metric=args.metric,
                       n_components=2)
        coords = reducer.fit(fps)

    else:
        raise ValueError(f"Unsupported method: {args.method}")

    with open(args.output, "w") as out:
        for mol_name, (x, y) in zip(mol_names, coords):
            x = str(round(x, 4))
            y = str(round(y, 4))
            out.write(mol_name + '\t' + x + '\t' + y + '\n')


if __name__ == "__main__":
    main()
