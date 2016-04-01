import mdtraj as md
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from fah_parameters import *

pdbid = "4IJ8"
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

pdb_filename = "./pdbs/%s_box.pdb" % pdbid

out_pdb_filename = "./equil/%s.pdb" % pdbid
final_step_pdb_filename = "./equil/%s_final_step.pdb" % pdbid
dcd_filename = "./equil/%s.dcd" % pdbid
log_filename = "./equil/%s.log" % pdbid

ff = app.ForceField(which_forcefield, which_water)

pdb = app.PDBFile(pdb_filename)

topology = pdb.topology
positions = pdb.positions

system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)

integrator = mm.LangevinIntegrator(temperature, equil_friction, equil_timestep)
system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)
print('Equilibrating...')

simulation.step(discard_steps)  # Don't even save the first XXX ps

simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.PDBReporter(out_pdb_filename, n_steps - 1))
simulation.reporters.append(app.StateDataReporter(open(log_filename, 'w'), output_frequency, step=True, time=True, speed=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(n_steps)

del simulation
del system
t = md.load(dcd_filename, top=out_pdb_filename)
t0 = t[-1]
t0.unitcell_lengths = t.unitcell_lengths.mean(0)
t0.save(out_pdb_filename)

del t
del t0
t = md.load(dcd_filename, top=out_pdb_filename)[-1]
t.save(final_step_pdb_filename)
