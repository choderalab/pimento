pimento
=======

##Manifest
* `Custom Chemistry` - contains custom chemistry definition used for covalently docking SGSmethyl inhibitor
* `SETD8 Crystal Structures` - original X-ray structures of SETD8 **[SETD8-DRG_refmac1 NOT YET PUBLISHED/CONFIDENTIAL]**
* ` SETD8 Crystal Structures` - structures docked in Schrodinger **[FILES LABELED NEW STRUCTURE CONFIDENTIAL]**
* `SETD8_newproject.prj` - Schrodinger Maestro directory
* `SETD8ligands_smiles` - positive hits and decoys provided by Luo lab **[NOT YET PUBLISHED/CONFIDENTIAL]**

##Schrodinger Directions

###Protein Prep Directions 
* In workspace, delete Chain B and any atoms labeled DU. 
* `Preprocess` using default settings + `Fill in missing side chains using Prime`
* On `Refine` tab, check `Sample water orientations` & `Use PROPKA ph=7.0` --> hit `optimize`
* On `Refine` tab, `minimize` using RMSD=.30A and `OPLS_2005` forcefield

###Ligand Preparation
* `Force Field` - OPLS_2005
* `Ionization` - Generate possible States at pH= 7.0 +/- 2.0 using Epik 
* Check `Desalt` & `Generate Tautomers`
* `Stereoisomers` - Used retain specified chiralities (vary other chiral centers) 
* `Generate low energy ring conformations` **1** per ligand`

###Covalent Docking Settings
* Ran 4 different docking job, ones for each ligand (SGSmethyl, SGSphenyl, Decoymethyl, Decoyphenyl)
* `Choose Reaction Type`: `Michael Addition` for *phenyl ligands and `Custom Chemistry` file from above for *methyl
* `Reactive Residue` is Cys311 
* `Docking box Center` Used SAM cofactor for `4IJ8` and the co-crystalized ligand in `newstructure`
* `Docking Mode` - Pose prediction (Thorough)
* `Output` - default (1 output pose, 1000 max number of ligands to report)



