#!/usr/bin/env python

# NOTE: Requires pdbfixer Python API and gaff2xml installed.

if __name__ == "__main__":

    import os
    import sys
    import glob
    import subprocess
    from argparse import ArgumentParser

    try:
        from gaff2xml import openeye
    except ImportError:
        print("Error: cannot import gaff2xml library!")
        sys.exit()

    try:
        from pdbfixer import PDBFixer
    except ImportError:
        print("Error: cannot import pdbfixer tool API!")
        sys.exit()

    from simtk.openmm.app import PDBFile

    # Constants
    #-----------

    # Input files
    RECEPTOR_FILENAME = 'receptor.pdb'
    SMILES_FILENAME = 'all_ligands.smi' # name of the input SMILES file
    LEAP_IN_FILENAME = 'leap.in' # tleap input file

    # Output files
    FIXED_RECEPTOR_FILENAME = 'receptor.pdbfixer.pdb' # name of the output file of PDBFixer
    TRIPOS_MOL2_FILENAME = 'ligand.tripos.mol2' # name of the output MOL2 file
    GAFF_MOL2_FILENAME = 'ligand.gaff.mol2' # name of the parametrized MOL2 file
    GAFF_FRCMOD_FILENAME = 'ligand.gaff.frcmod'
    LEAP_OUT_FILENAME = 'leap.out' # tleap output file

    # Argument parsing
    #------------------

    parser = ArgumentParser(description='Convert the SMILES file into the input file need by Yank'
        ' and parametrize the molecule.')
    parser.add_argument('-nmol', default=0, type=int, help='The line of the SMILES file describing the molecule'
        ' to analyze.', dest='n_mol')
    args = parser.parse_args()

    # SMILES to MOL2 conversion
    #---------------------------

    print('Generating ligand Tripos mol2 file from SMILES...')

    # Converting n_mol-th string of input file from SMILES format to OEMol
    smiles_file = open(SMILES_FILENAME, 'r')
    mol = openeye.smiles_to_oemol(smiles_file.readlines()[args.n_mol].strip())
    smiles_file.close()

    # Assigning charges to molecule
    charged_mol = openeye.get_charges(mol)

    # Converting OEMol to MOL2 format
    openeye.molecule_to_mol2(charged_mol, tripos_mol2_filename=TRIPOS_MOL2_FILENAME)

    # Parametrization
    #-----------------

    print('Parametrizing ligand with GAFF and AM1-BCC carges...')

    bash_cmd = 'antechamber'
    bash_cmd += ' -fi mol2 -i ' + TRIPOS_MOL2_FILENAME # antechamber input file in mol2 format
    bash_cmd += ' -fo mol2 -o ' + GAFF_MOL2_FILENAME # antechamber output file in mol2 format
    subprocess.call(bash_cmd.split())

    bash_cmd = 'parmchk'
    bash_cmd += ' -f mol2 -i ' + GAFF_MOL2_FILENAME #parmchk input file in mol2 format
    bash_cmd += ' -o ' + GAFF_FRCMOD_FILENAME # parmchk output file in frcmod format
    subprocess.call(bash_cmd.split())

    # Fixing receptor structure
    #---------------------------

    print("Preparing receptor by adding missing atoms...")

    fixer = PDBFixer(filename=RECEPTOR_FILENAME)

    # Fill in missing atoms
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Expunge waters and ions
    fixer.removeHeterogens(False)

    # Write pdb file
    PDBFile.writeFile(fixer.topology, fixer.positions, open(FIXED_RECEPTOR_FILENAME, 'w'))

    # AMBER prmtop/inpcrd files
    #---------------------------

    print("Creating AMBER prmtop/inpcrd files...")

    bash_cmd = 'tleap -f ' + LEAP_IN_FILENAME
    leap_out_file = open(LEAP_OUT_FILENAME, 'w')
    subprocess.call(bash_cmd.split(), stdout=leap_out_file)
    leap_out_file.close()
