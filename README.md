### RDKit scripts
###### let's make life easier

The purpose of this repository is to collect useful scripts which mainly use RDKit. Contributions are welcome!

#### Comments and recommendations for contributors:
1. There is a `read_input.py` script which contains the function `read_input`. It reads molecules from SMI, SDF, SDF.GZ and PKL (pickled molecules as tuples of mol and mol_title) files and STDIN (SMI and SDF formats are supported) and it returns tuples of (mol, mol_title). This is a generator and can be applied to process large collections of molecules. I advise to use this function if you do not need other data from input files.
2. There is `_template.py` file which can be used as a template for new scripts. Please do not change names for input, output, ncpu and verbose arguments. This will help to make command line arguments consistent across scripts.
3. Add help messages to your scripts.
4. Ideally scripts should be able to communicate with STDIN and STDOUT to combine them with pipes. I implemented this in `gen_stereo_rdkit.py` and `gen_conf_rdkit.py`.
5. All scripts can contain errors, so use them on your own risk. If you will find a mistake please create the issue and we will fix it. However, we constantly revise old scripts and fix errors because every found mistake is penultimate.

#### Particular scripts

`vina_dock.py` - performs parallel docking in fully automatic way: protonation of ligands, generation of 3D PDBQT representations, docking and storage of output in a single DB. It can run on a single machine as well as on distributed infrastructure using dask, interrupted computations can be easily resumed. It has several dependencies. See help for more details and examples.

##### Happy RDKiiting! :)
