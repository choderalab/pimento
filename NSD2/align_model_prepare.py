## Rafal P. Wiewiora, ChoderaLab, @rafwiewiora March 2016
from simtk.openmm import app
import modeller
import modeller.automodel
import os
import glob
from fah_parameters import *

## Note this was developed first for OpenMM-7.0beta + you need the ZNII-multiSite.xml ffxml
# an updated version of that ffxml will be distributed with OpenMM as standard starting 7.1


if not os.path.exists('pdb_outputs/'):
    os.mkdir('pdb_outputs')
if not os.path.exists('modeller_outputs/'):
    os.mkdir('modeller_outputs')

pdb_inputs = ('3ooi.pdb', '4yz8.pdb')
pdb_input_codes = ('3ooi', '4yz8')
target_pir = 'input_files/nsd2.pir'
target_code = 'NSD2'

# remove all but protein
for pdb_name in pdb_inputs:
    pdb = app.PDBFile('input_files/' + pdb_name)
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.delete([res for res in mod.topology._chains[1].residues()])
    pdb_out_name = 'pdb_outputs/' + pdb_name
    app.PDBFile.writeFile(mod.topology, mod.positions, open(pdb_out_name, 'w'))

# make alignments
env = modeller.environ()
aln = modeller.alignment(env)
models = dict()
for pdb_code in pdb_input_codes:
    models[pdb_code] = modeller.model(env, file=('pdb_outputs/%s.pdb' % pdb_code),
                             model_segment=('FIRST:A','LAST:A'))
    aln.append_model(models[pdb_code], align_codes=pdb_code,
                     atom_files=('pdb_outputs/%s.pdb' % pdb_code))
aln.append(file=target_pir, align_codes=target_code)
aln.align2d()
aln.write(file=('modeller_outputs/%s.ali' % target_code), alignment_format='PIR')

#add chain break and 3 dots to the aligns for 3X ZN
ali_file = open('modeller_outputs/%s.ali' % target_code)
lines = [line for line in ali_file]
ali_file.close()
new_lines = []
for line in lines:
    # * end of sequence, add /... before
    if '*' in line:
        star_index = line.index('*')
        new_line = line[:star_index] + '/...' + '*\n'
    # + starts 3-digit (hardcoded here) number of AAs - increase by 3
    elif '+' in line:
        plus_index = line.index('+')
        aa_number = int(line[plus_index+1:plus_index+4])
        aa_number += 3
        new_line = line[:plus_index] + '+' + str(aa_number) + line[plus_index+4:]
    else:
        new_line = line
    new_lines.append(new_line)
f = open('modeller_outputs/%s_withZN.ali' % target_code, 'w')
for write_line in new_lines:
    f.write(write_line)
f.close()

# here we regenrate the pdb_outputs files that are input to modelling - we want
# ZN's in for modelling now, but didn't for aligmment - we added to dots separately
for pdb_name in pdb_inputs:
    pdb = app.PDBFile('input_files/' + pdb_name)
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.delete([res for res in mod.topology._chains[1].residues()
                if res.name != 'ZN'])
    # have ZN's and protein in the same chain to avoid problems with alignment
    for res in mod.topology._chains[1].residues():
        mod.topology._chains[0]._residues.append(res)
    mod.topology._chains = [mod.topology._chains[0]]
    pdb_out_name = 'pdb_outputs/' + pdb_name
    app.PDBFile.writeFile(mod.topology, mod.positions, open(pdb_out_name, 'w'))

# do modelling
env = modeller.environ()
# must turn hetatm to True for ZN's
env.io.hetatm = True
a = modeller.automodel.automodel(env,
                       alnfile=('modeller_outputs/%s_withZN.ali' % target_code),
                       knowns=pdb_input_codes, sequence=target_code)
a.starting_model= 1
a.ending_model  = 1
a.make()

# move the outputs to modeller_outputs
modelling_outputs = glob.glob('%s*' % target_code)
for output in modelling_outputs:
    os.system('mv %s modeller_outputs/' % output)

# do app.modeller stuff - hydrogens (pay attention to deprotonation), water, EP
print('Loading into OpenMM for final prep.')
pdb = app.PDBFile('modeller_outputs/%s.B99990001.pdb' % target_code)
mod = app.Modeller(pdb.topology, pdb.positions)
# indexes of cysteines that need deprotonation
# IDENTIFICATION OF CYSTEINES DONE VISUALLY BY RPW IN PYMOL
# TODO: DO THE CYSTEINE IDENTIFICATIONS BY ALGORITHM (DISTANCE)
cysteines = (172, 219, 221, 226, 74, 80, 69, 60, 54, 44, 46) # these are 1-indexed!!
ff = app.ForceField('amber99sbildn.xml', 'ZNII-multiSite.xml', 'tip3p.xml')
variants = [None] * len(list(mod.topology.residues()))
for cys in cysteines:
    variants[cys-1] = 'CYX' # remember CYM is the same as CYX at this stage - no hydrogen

print('Adding hydrogens.')
mod.addHydrogens(variants=variants)
print('Adding extra particles')
mod.addExtraParticles(ff)
print('Adding solvent')
mod.addSolvent(ff, padding=padding)
#write back to pdb for taking to the cluster
print('Saving PDB.')
app.PDBFile.writeFile(mod.topology, mod.positions, open('pdb_outputs/%s.pdb' % target_code, 'w'))
