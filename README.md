### RDKit scripts
###### let's make life easier

The purpose of this repository is to collect useful scripts which mainly use RDKit. Contributions are welcome!

Some scripts may require further dependencies.

#### Comments and recommendations for contributors:
1. There is a `read_input.py` script which contains the function `read_input`. It reads molecules from SMI, SDF, SDF.GZ and PKL (pickled molecules as tuples of mol and mol_title) files and STDIN (SMI and SDF formats are supported) and it returns tuples of (mol, mol_title). This is a generator and can be applied to process large collections of molecules. I advise to use this function if you do not need other data from input files.
2. There is `_template.py` file which can be used as a template for new scripts. Please do not change names for input, output, ncpu and verbose arguments. This will help to make command line arguments consistent across scripts.
3. Add help messages to your scripts.
4. Ideally scripts should be able to communicate with STDIN and STDOUT to combine them with pipes. I implemented this in `gen_stereo_rdkit.py` and `gen_conf_rdkit.py`.
5. All scripts can contain errors, so use them on your own risk. If you will find a mistake please create the issue and we will fix it. However, we constantly revise old scripts and fix errors because every found mistake is penultimate.

#### Particular scripts

##### Manipulate with SDF:  
  
`add_prefix` - add a prefix to molecule names in SDF file.  

`extractsdf` - extract molecule names and field values from input SDF.  
  
`extract_mol_by_name` - extract molecules by name (partial name matching) to new SDF file   

`insert_sdf` - add data from a text file as additional fields to input SDF file  
  
`remove_dupl_by_field` - remove entries from SDF file having duplicated mol title or field value.  
  
`rename_mols` - identify identical entries in SDF (conformers) and rename them in identical manner. 

`sdf_field2title` - insert values of a given SDF field into molecular title.  
  
`strip_blank_lines` - remove empty lines in multi-line field values in input SDF.  
  
##### Format and file (inter)conversion:  

`cansmi` - return canonical SMILES of input molecules.  
  
`frags2mols` - save disconnected components of input molecules as individual molecules with added suffix to the name.  
  
`molchemaxon2pdb` - convert molecules from the input file to separate pdb files. Conformer generation is performed by RDKit. Major tautomeric form at the given pH is generated by ChemAxon.  
  
`mols2pdb` - convert input molecules from SMILES or SDF file to individual PDB files. Hydrogens will be added and a random conformer will be generated if the molecule does not have 3D coordinates.  
  
`pkl2sdf` - convert PKL file to SDF file. Specifically useful for conversion of generated conformers stored in PKL format by `gen_conf_rdkit`.  
  
`sdf2mols` - split SDF to multiple MOL files.  
  
`sdf2pkl` - convert SDF to multi-conformer PKL file. Conformers are recognized by mol title and should go sequentially in input SDF.  
  
`smi2sdf` - Convert SMILES to SDF including additional fields if they are named and exist in SMILES file  
  
`split_pdb` - split PDB by chains and save to separate PDB files.  

##### Manipulate with Mol objects (calc properties, generate conformers/stereoisomers, filter compounds, etc):  

`add_h` - hydrogenize input files  
  
`calc_center_rdkit` - calculate the center of coordinates of all atoms in a molecule(s).  
  
`count_undefined_stereocenters` - return to STDOUT names and the number of undefined stereocenters in input molecules.  
  
`discard_compounds_rdkit` - remove multi-component compounds and compounds with non-organic atoms.  
  
`draw_mols` - return PNG images of input molecules.  
  
`filter_conf` - filter conformers by RMS value.  
  
`gen_conf_rdkit` - generate conformers.  
  
`gen_stereo_rdkit` - enumerate stereoisomers (tetrahedral and double bond).  
  
`get_mol_center.py` - returns a geometrical center of a molecule
  
`get_substr` - filter input molecules by SMARTS, multiple SMARTS are allowed, negative matching is possible.  
  
`get_total_charge` - calculate total formal charge of input MOL files.  
  
`keep_largest` - keep the largest fragment by the number of heavy atoms in each compound record. If components have the same number a random one will be selected.  
  
`mirror_mols.py` - return mirrored 3D input structure and optionally rename it. Useful for generation of enantiomers of molecules with axial/planar chirality.
  
`murcko` - return Murcko scaffolds ignoring stereoconfiguration.  
  
`physchem_calc` - calculate various physicochemical properties of input molecules (MW, logP, TPSA, QED, etc).  
  
`pmapper_descriptors` - calculate 3D pharmacophore descriptors with `pmapper` and remove rarely occurred ones. Useful for QSAR modeling.  
  
`remove_dupl_rdkit` - remove duplicates by InChi keys comparison within the input file or relatively to a reference file.  
  
`rmsd_rdkit` - calculate RMSD for input MOL2/PDBQT/SDF files. Automatically calculate RMSD for maximum common substructure if full atom matching was failed. Symmetry checking was implemented.  
  
`sanitize_rdkit` - remove compounds with RDKit sanitization errors and add to output molecules the number of double bonds, unspecified stereocenters and total charge.  

`sphere_exclusion` - return names of a diverse subset of input compounds  
  
`test_pains` - return a list of SMILES matched PAINS.  

##### Supplementary scripts:

`binning` - take a table with variable values and return a table with binned values according to supplied thresholds.  

##### Removed scripts
`vina_dock` has been removed because there is a separate repository for automatic and distributed docking which includes docking with Vina, smina and gnina - https://github.com/ci-lab-cz/easydock.

##### Happy RDKiiting! :)
