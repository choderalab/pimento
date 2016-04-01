import os
import numpy as np
from fah_parameters import *
import simtk.openmm.app as app
import pdbfixer
import mdtraj as md

PADDING = padding  # From fah_parameters.py

def sort_atoms(t0):
    """Use lexical sort to standardize atom order, returning a new trajectory."""
    top0, bonds0 = t0.top.to_dataframe()

    top0["x"] = t0.xyz[0, :, 0]
    top0["y"] = t0.xyz[0, :, 1]
    top0["z"] = t0.xyz[0, :, 2]

    top0 = top0.sort_index(by=["resSeq", "name"])

    n_atoms = len(top0)
    top0.serial = np.arange(1, n_atoms + 1)

    top0.index = np.arange(n_atoms)

    bonds = np.zeros((0, 2), 'int')
    top = md.Topology.from_dataframe(top0, bonds)
    xyz = np.array([top0.x.values, top0.y.values, top0.z.values]).T[None]
    traj = md.Trajectory(xyz, top)

    return traj



def fix(pdbid, missing_residues, padding=PADDING, mutations=None):
    fixer = pdbfixer.PDBFixer(filename="%s.pdb" % pdbid)

    if mutations is not None:
        fixer.applyMutations(mutations[0], mutations[1])

    fixer.missingResidues = missing_residues
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    numChains = len(list(fixer.topology.chains()))
    fixer.removeChains(range(1, numChains))
    file_handle = open("./pdbs/%s_fixed0.pdb" % pdbid, 'wb')
    app.PDBFile.writeFile(fixer.topology, fixer.positions, file_handle)
    file_handle.close()
    
    traj0 = md.load("./pdbs/%s_fixed0.pdb" % pdbid)
    traj = sort_atoms(traj0)
    filename1 = "./pdbs/%s_fixed1.pdb" % pdbid
    traj.save(filename1)
    
    filename = "./pdbs/%s_fixed.pdb" % pdbid
    # Need to sync the protonation state of one histidine
    
    cmd = """grep  -v 'HE1 HIS A 119' %s |grep  -v 'HE2 HIS A 119'|grep  -v 'HD1 HIS A 119'|grep  -v 'HD2 HIS A 119' > %s""" % (filename1, filename)
    os.system(cmd)
    #os.system("grep  -v 'HE2 HIS A 119' %s > %s""" % (filename1, filename))
    #os.system("grep  -v 'HD1 HIS A 119' %s > %s""" % (filename1, filename))
    #os.system("grep  -v 'HD2 HIS A 119' %s > %s""" % (filename1, filename))

    ff_name = "amber99sbildn"
    water_name = 'tip3p'

    which_forcefield = "%s.xml" % ff_name
    which_water = '%s.xml' % water_name

    out_pdb_filename = "./equil/equil.pdb"
    ff = app.ForceField(which_forcefield, which_water)
    
    pdb = app.PDBFile(filename)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    variants = [None for i in range(161)]
    variants[118] = "HIE"
    
    modeller.addHydrogens(ff, variants=variants)
    modeller.addSolvent(ff, padding=padding)

    app.PDBFile.writeFile(modeller.topology, modeller.positions, open("./pdbs/%s_box.pdb" % pdbid, 'w'))

fix("1ZKK", {(0, 0):["ARG"]})  # Add the N terminal ARG for 1ZKK
mutations=(["SER-343-CYS"], "A")
#mutations = None
fix("4IJ8", {}, mutations=mutations)  # Don't add any residues for 4IJ8

# NOTE: this relies on the fact that waters are placed last.

pdbids = ["1ZKK", "4IJ8"]
trajectories = [md.load("././pdbs/%s_box.pdb" % pdbid) for pdbid in pdbids]
n_atoms = min([t.n_atoms for t in trajectories])

for k, t in enumerate(trajectories):
    pdbid = pdbids[k]
    t.atom_slice(np.arange(n_atoms)).save("./pdbs/%s_box.pdb" % pdbid)
