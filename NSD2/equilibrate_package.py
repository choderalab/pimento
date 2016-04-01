from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from fah_parameters import *
import mdtraj as md

pdbid = 'NSD2'

pdb_filename = "./pdb_outputs/%s.pdb" % pdbid
out_pdb_filename = "./equil/%s.pdb" % pdbid
dcd_filename = "./equil/%s.dcd" % pdbid
log_filename1 = "./equil/%s.log" % (pdbid + '_1fs')
log_filename2 = "./equil/%s.log" % (pdbid + '_2fs')

system_filename = "./package/system.xml"
integrator_filename = "./package/integrator.xml"
state_filename = "./package/state.xml"

print("Loading PDB...")
pdb = PDBFile(pdb_filename, extraParticleIdentifier='') # extraParticleIdentifier for loading dummy atoms

topology = pdb.topology
positions = pdb.positions
        
ff = ForceField('amber99sbildn.xml', 'ZNII-multiSite.xml', 'tip3p.xml')       
platform = Platform.getPlatformByName('CUDA')

# need to minimize with CutoffPeriodic nonbondedMethod to avoid blowing up
print("Preparing minimization system...")        
system = ff.createSystem(topology, nonbondedMethod=CutoffPeriodic, nonbondedCutoff=cutoff, constraints=HBonds)        

print("Preparing simulation for minimization...")
integrator = VerletIntegrator(equil_timestep)
simulation = app.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)
print("Initial energy is %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print("Energy after minimization is %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
positions = simulation.context.getState(getPositions=True).getPositions()
vectors = simulation.context.getState().getPeriodicBoxVectors()
del simulation

print("Preparing new 1fs system...")
system = ff.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=cutoff, constraints=HBonds)
system.setDefaultPeriodicBoxVectors(vectors[0], vectors[1], vectors[2])        
integrator = LangevinIntegrator(temperature, equil_friction, equil_timestep)
system.addForce(MonteCarloBarostat(pressure, temperature, barostat_frequency))

print("Preparing 1fs simulation...")
simulation = app.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

print("Initial energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print("Energy after minimization is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))

simulation.context.setVelocitiesToTemperature(temperature)

print("Running simulation - discard_steps...")
simulation.step(discard_steps)

print("Appending StateDataReporter...")
simulation.reporters.append(app.StateDataReporter(open(log_filename1, 'w'), output_frequency, step=True, time=True, speed=True, potentialEnergy=True, temperature=True, density=True))

print("Running production steps with 1fs timestep...")
simulation.step(n_steps)   
print('Done with 1fs equilibration!')
print("Final energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
positions = simulation.context.getState(getPositions=True).getPositions()
vectors = simulation.context.getState().getPeriodicBoxVectors()
del simulation

print("Preparing new integrator and updating system with new PeriodicBoxVectors for 2fs equilibration...")
integrator = LangevinIntegrator(temperature, equil_friction, timestep)
system.setDefaultPeriodicBoxVectors(vectors[0], vectors[1], vectors[2]) 

print("Preparing 2fs simulation...")
simulation = app.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)

print("Initial energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print("Energy after minimization is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))

simulation.context.setVelocitiesToTemperature(temperature)

print("Appending reporters...")
simulation.reporters.append(app.DCDReporter(dcd_filename, output_frequency))
simulation.reporters.append(app.PDBReporter(out_pdb_filename, n_steps - 1))
simulation.reporters.append(app.StateDataReporter(open(log_filename2, 'w'), output_frequency, step=True, time=True, speed=True, potentialEnergy=True, temperature=True, density=True))

print("Running production steps with 2fs timestep...")
simulation.step(n_steps)   
print('Done with 2fs equilibration!')
print("Final energy is: %s" % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
positions = simulation.context.getState(getPositions=True).getPositions()
vectors = simulation.context.getState().getPeriodicBoxVectors()
del simulation

# Do packaging
print('Initiating packaging...')
def write_file(filename, contents):
    with open(filename, 'w') as outfile:
        outfile.write(contents)
        
print('Preparing new integrator and updating system with new PeriodicBoxVectors...')
system.setDefaultPeriodicBoxVectors(vectors[0], vectors[1], vectors[2]) 
integrator = LangevinIntegrator(temperature, friction, timestep)

print("Preparing new simulation...")
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(output_frequency)

print("Writing system and integrator...")
write_file(system_filename, XmlSerializer.serialize(system))
write_file(integrator_filename, XmlSerializer.serialize(integrator))

print("Writing state...")
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
write_file(state_filename, XmlSerializer.serialize(state))
print('All done!')
        
        
        
        
