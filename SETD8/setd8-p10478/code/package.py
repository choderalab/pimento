import os
import time
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u
from fah_parameters import *

def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)

pdbids = ["1ZKK", "4IJ8"]
ff_name = "amber99sbildn"
water_name = 'tip3p'

which_forcefield = "%s.xml" % ff_name
which_water = '%s.xml' % water_name

for run, pdbid in enumerate(pdbids):
    run_dir =  os.path.join(os.path.expanduser("~"), "dat/fah_data/input_data/10478/RUNS/RUN%d/" % run)

    nclones = 500

    system_filename = os.path.join(run_dir, "system.xml")
    integrator_filename = os.path.join(run_dir, "integrator.xml")

    pdb_filename = "./equil/%s_final_step.pdb" % pdbid

    pdb = app.PDBFile(pdb_filename)
    topology = pdb.topology
    positions = pdb.positions

    ff = app.ForceField(which_forcefield, which_water)
    system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=app.HBonds)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, barostat_frequency))

    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature)

    write_file(system_filename, mm.XmlSerializer.serialize(system))
    write_file(integrator_filename, mm.XmlSerializer.serialize(integrator))

    for clone_index in range(nclones):
        print(run, clone_index)
        simulation.context.setVelocitiesToTemperature(temperature)
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
        state_filename = os.path.join(run_dir, 'state%d.xml' % clone_index)
        serialized = mm.XmlSerializer.serialize(state)
        write_file(state_filename, serialized)
