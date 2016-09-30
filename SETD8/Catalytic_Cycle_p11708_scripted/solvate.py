from simtk.openmm.app import PDBFile, Modeller, ForceField
from simtk import unit as u

ff_sam = ForceField('amber99sbildn.xml', 'tip3p.xml', 'parameters/gaff.xml', 'parameters/ffptm.xml', 'parameters/SAM.xml') 
ff_sah = ForceField('amber99sbildn.xml', 'tip3p.xml', 'parameters/gaff.xml', 'parameters/ffptm.xml', 'parameters/SAH.xml') 

# load all 9 structurres in
set8 = PDBFile('no_solvent/SET8.pdb')
set8_p = PDBFile('no_solvent/SET8_P.pdb')
set8_mep = PDBFile('no_solvent/SET8_MeP.pdb')

set8_sam = PDBFile('no_solvent/SET8_SAM.pdb')
set8_sah = PDBFile('no_solvent/SET8_SAH.pdb')

set8_p_sam = PDBFile('no_solvent/SET8_P_SAM.pdb')
set8_p_sah = PDBFile('no_solvent/SET8_P_SAH.pdb')

set8_mep_sam = PDBFile('no_solvent/SET8_MeP_SAM.pdb')
set8_mep_sah = PDBFile('no_solvent/SET8_MeP_SAH.pdb')

# create Modeller objects
set8_modeller = Modeller(set8.topology, set8.positions)
set8_p_modeller = Modeller(set8_p.topology, set8_p.positions)
set8_mep_modeller = Modeller(set8_mep.topology, set8_mep.positions)

set8_sam_modeller = Modeller(set8_sam.topology, set8_sam.positions)
set8_sah_modeller = Modeller(set8_sah.topology, set8_sah.positions)

set8_p_sam_modeller = Modeller(set8_p_sam.topology, set8_p_sam.positions)
set8_p_sah_modeller = Modeller(set8_p_sah.topology, set8_p_sah.positions)

set8_mep_sam_modeller = Modeller(set8_mep_sam.topology, set8_mep_sam.positions)
set8_mep_sah_modeller = Modeller(set8_mep_sah.topology, set8_mep_sah.positions)

# solvate
set8_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_p_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_mep_modeller.addSolvent(ff_sam, padding=1*u.nanometers)

set8_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

set8_p_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_p_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

set8_mep_sam_modeller.addSolvent(ff_sam, padding=1*u.nanometers)
set8_mep_sah_modeller.addSolvent(ff_sah, padding=1*u.nanometers)

# save all in with_solvent
PDBFile.writeFile(set8_modeller.topology, set8_modeller.positions, open('with_solvent/SET8.pdb', 'w'))
PDBFile.writeFile(set8_p_modeller.topology, set8_p_modeller.positions, open('with_solvent/SET8_P.pdb','w'))
PDBFile.writeFile(set8_mep_modeller.topology, set8_mep_modeller.positions, open('with_solvent/SET8_MeP.pdb', 'w'))

PDBFile.writeFile(set8_sam_modeller.topology, set8_sam_modeller.positions, open('with_solvent/SET8_SAM.pdb','w'))
PDBFile.writeFile(set8_sah_modeller.topology, set8_sah_modeller.positions, open('with_solvent/SET8_SAH.pdb', 'w'))

PDBFile.writeFile(set8_p_sam_modeller.topology, set8_p_sam_modeller.positions, open('with_solvent/SET8_P_SAM.pdb','w'))
PDBFile.writeFile(set8_p_sah_modeller.topology, set8_p_sah_modeller.positions, open('with_solvent/SET8_P_SAH.pdb', 'w'))

PDBFile.writeFile(set8_mep_sam_modeller.topology, set8_mep_sam_modeller.positions, open('with_solvent/SET8_MeP_SAM.pdb','w'))
PDBFile.writeFile(set8_mep_sah_modeller.topology, set8_mep_sah_modeller.positions, open('with_solvent/SET8_MeP_SAH.pdb', 'w'))
