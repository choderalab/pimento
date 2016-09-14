import mdtraj as md
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

trajs = glob.glob('../data/*.h5')

# load structures and get atom indices for appropriate selections
data = md.load('../data/run0-clone0.h5')[0]
zkk = md.load_pdb('files/pdb1zkk.ent')
ij8 = md.load_pdb('files/pdb4ij8.ent')
apo = md.load('files/Apo.pdb')
inhibitor = md.load('files/Inhibitor.pdb')

data_selection = data.top.select('not (resid 0 or resid 159 or resid 160 or resid 76) and name CA')
zkk_selection = zkk.top.select('chainid 0 and not (resid 158 or resid 159 or resid 75) and name CA')
ij8_selection = ij8.top.select('chainid 0 and not (resid 0 or resid 159 or resid 160 or resid 76) and name CA')
apo_selection = apo.top.select('chainid 0 and not (resid 0 or resid 1 or resid 2 or resid 3) and name CA')
inhibitor_selection = inhibitor.top.select('chainid 0 and not (resid 158 or resid 159 or resid 75) and name CA')

zkk_ca = zkk.atom_slice(zkk_selection)
ij8_ca = ij8.atom_slice(ij8_selection)
apo_ca = apo.atom_slice(apo_selection)
inhibitor_ca = inhibitor.atom_slice(inhibitor_selection)

# calculate rmsds - save .npy's
zkk_rmsds = []
ij8_rmsds = []
apo_rmsds = []
inhibitor_rmsds = []

current_trajectory = 1

for traj in trajs:
    print('Currently calculating for trajectory: %d' % current_trajectory)
    current_trajectory += 1
    traj_ = md.load(traj)
    traj_ca = traj_.atom_slice(data_selection)
    zkk_rmsd = md.rmsd(traj_ca, zkk_ca)
    zkk_rmsds.append(zkk_rmsd)
    ij8_rmsd = md.rmsd(traj_ca, ij8_ca)
    ij8_rmsds.append(ij8_rmsd)
    apo_rmsd = md.rmsd(traj_ca, apo_ca)
    apo_rmsds.append(apo_rmsd)
    inhibitor_rmsd = md.rmsd(traj_ca, inhibitor_ca)
    inhibitor_rmsds.append(inhibitor_rmsd)

# save results
np.save('CA_rmsds/1zkk_rmsds.npy', zkk_rmsds)
np.save('CA_rmsds/4ij8_rmsds.npy', ij8_rmsds)
np.save('CA_rmsds/apo_rmsds.npy', apo_rmsds)
np.save('CA_rmsds/inhibitor_rmsds.npy', inhibitor_rmsds)

# make plots
plt.figure()
for i in range(len(zkk_rmsds)):
    plt.plot(zkk_rmsds[i])
plt.title('CA RMSDs to 1ZKK - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_rmsds.pdf')

plt.figure()
for i in range(len(ij8_rmsds)):
    plt.plot(ij8_rmsds[i])
plt.title('CA RMSDs to 4IJ8 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_rmsds.pdf')

plt.figure()
for i in range(len(apo_rmsds)):
    plt.plot(apo_rmsds[i])
plt.title('CA RMSDs to Apo - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')  
plt.savefig('CA_rmsds/apo_rmsds.pdf')

plt.figure()
for i in range(len(inhibitor_rmsds)):
    plt.plot(inhibitor_rmsds[i])
plt.title('CA RMSDs to Inhibitor - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')  
plt.savefig('CA_rmsds/inhibitor_rmsds.pdf')

# now do the calculations BY DOMAIN
# ranges - see README for more info
# N-flanking: 0 - 22 (2 first missing)
# SET(1): 23 - 56
# SET-I: 57 - 92
# SET(2): 93 - 142
# C-flanking: 143 - 157 (1 last missing)

zkk_Nflank_rmsds = []
ij8_Nflank_rmsds = []
apo_Nflank_rmsds = []
inhibitor_Nflank_rmsds = []

zkk_SET1_rmsds = []
ij8_SET1_rmsds = []
apo_SET1_rmsds = []
inhibitor_SET1_rmsds = []

zkk_SETI_rmsds = []
ij8_SETI_rmsds = []
apo_SETI_rmsds = []
inhibitor_SETI_rmsds = []

zkk_SET2_rmsds = []
ij8_SET2_rmsds = []
apo_SET2_rmsds = []
inhibitor_SET2_rmsds = []

zkk_Cflank_rmsds = []
ij8_Cflank_rmsds = []
apo_Cflank_rmsds = []
inhibitor_Cflank_rmsds = []

# zkk starts at resid 0, ij8 at 1 - add 1, apo at 4 - add 4, inhibitor at 0
zkk_ca_Nflank = zkk_ca.atom_slice(zkk_ca.top.select('resid >= 0 and resid <= 22'))
ij8_ca_Nflank = ij8_ca.atom_slice(ij8_ca.top.select('resid >= 1 and resid <= 23'))
apo_ca_Nflank = apo_ca.atom_slice(apo_ca.top.select('resid >= 4 and resid <= 26'))
inhibitor_ca_Nflank = inhibitor_ca.atom_slice(inhibitor_ca.top.select('resid >= 0 and resid <= 22'))

zkk_ca_SET1 = zkk_ca.atom_slice(zkk_ca.top.select('resid >= 23 or resid <= 56'))
ij8_ca_SET1 = ij8_ca.atom_slice(ij8_ca.top.select('resid >= 24 or resid <= 57'))
apo_ca_SET1 = apo_ca.atom_slice(apo_ca.top.select('resid >= 27 or resid <= 60'))
inhibitor_ca_SET1 = inhibitor_ca.atom_slice(inhibitor_ca.top.select('resid >= 23 or resid <= 56'))

zkk_ca_SETI = zkk_ca.atom_slice(zkk_ca.top.select('resid >= 57 or resid <= 92'))
ij8_ca_SETI = ij8_ca.atom_slice(ij8_ca.top.select('resid >= 58 or resid <= 93'))
#did 96 to 95
apo_ca_SETI = apo_ca.atom_slice(apo_ca.top.select('resid >= 61 or resid <= 95'))
inhibitor_ca_SETI = inhibitor_ca.atom_slice(inhibitor_ca.top.select('resid >= 57 or resid <= 92'))

zkk_ca_SET2 = zkk_ca.atom_slice(zkk_ca.top.select('resid >= 93 or resid <= 142'))
ij8_ca_SET2 = ij8_ca.atom_slice(ij8_ca.top.select('resid >= 94 or resid <= 143'))
# did 97 to 96 and 146 to 145
apo_ca_SET2 = apo_ca.atom_slice(apo_ca.top.select('resid >= 96 or resid <= 145'))
inhibitor_ca_SET2 = inhibitor_ca.atom_slice(inhibitor_ca.top.select('resid >= 93 or resid <= 142'))

zkk_ca_Cflank = zkk_ca.atom_slice(zkk_ca.top.select('resid >= 143'))
ij8_ca_Cflank = ij8_ca.atom_slice(ij8_ca.top.select('resid >= 144'))
# did 147 to 146
apo_ca_Cflank = apo_ca.atom_slice(apo_ca.top.select('resid >= 146'))
inhibitor_ca_Cflank = inhibitor_ca.atom_slice(inhibitor_ca.top.select('resid >= 143'))

current_trajectory = 1

for traj in trajs:
    print('Currently calculating for trajectory: %d' % current_trajectory)
    current_trajectory += 1
    traj_ = md.load(traj)
    
    traj_ca = traj_.atom_slice(data_selection)
    # add 1 to resid - they stay same from original traj and we started at 1
    traj_ca_Nflank = traj_ca.atom_slice(traj_ca.top.select('resid >= 1 and resid <= 23'))
    traj_ca_SET1 = traj_ca.atom_slice(traj_ca.top.select('resid >= 24 or resid <= 57'))
    traj_ca_SETI = traj_ca.atom_slice(traj_ca.top.select('resid >= 58 or resid <= 93'))
    traj_ca_SET2 = traj_ca.atom_slice(traj_ca.top.select('resid >= 94 or resid <= 143'))
    traj_ca_Cflank = traj_ca.atom_slice(traj_ca.top.select('resid >= 144'))
    
    #N_flank rmsds
    zkk_rmsd = md.rmsd(traj_ca_Nflank, zkk_ca_Nflank)
    ij8_rmsd = md.rmsd(traj_ca_Nflank, ij8_ca_Nflank)
    apo_rmsd = md.rmsd(traj_ca_Nflank, apo_ca_Nflank)
    inhibitor_rmsd = md.rmsd(traj_ca_Nflank, inhibitor_ca_Nflank)
    zkk_Nflank_rmsds.append(zkk_rmsd)
    ij8_Nflank_rmsds.append(ij8_rmsd)
    apo_Nflank_rmsds.append(apo_rmsd)
    inhibitor_Nflank_rmsds.append(inhibitor_rmsd)
    
    # SET1
    zkk_rmsd = md.rmsd(traj_ca_SET1, zkk_ca_SET1)
    ij8_rmsd = md.rmsd(traj_ca_SET1, ij8_ca_SET1)
    apo_rmsd = md.rmsd(traj_ca_SET1, apo_ca_SET1)
    inhibitor_rmsd = md.rmsd(traj_ca_SET1, inhibitor_ca_SET1)
    zkk_SET1_rmsds.append(zkk_rmsd)
    ij8_SET1_rmsds.append(ij8_rmsd)
    apo_SET1_rmsds.append(apo_rmsd)
    inhibitor_SET1_rmsds.append(inhibitor_rmsd)
    
    # SETI
    zkk_rmsd = md.rmsd(traj_ca_SETI, zkk_ca_SETI)
    ij8_rmsd = md.rmsd(traj_ca_SETI, ij8_ca_SETI)
    apo_rmsd = md.rmsd(traj_ca_SETI, apo_ca_SETI)
    inhibitor_rmsd = md.rmsd(traj_ca_SETI, inhibitor_ca_SETI)
    zkk_SETI_rmsds.append(zkk_rmsd)
    ij8_SETI_rmsds.append(ij8_rmsd)
    apo_SETI_rmsds.append(apo_rmsd)
    inhibitor_SETI_rmsds.append(inhibitor_rmsd)
    
    # SET2
    zkk_rmsd = md.rmsd(traj_ca_SET2, zkk_ca_SET2)
    ij8_rmsd = md.rmsd(traj_ca_SET2, ij8_ca_SET2)
    apo_rmsd = md.rmsd(traj_ca_SET2, apo_ca_SET2)
    inhibitor_rmsd = md.rmsd(traj_ca_SET2, inhibitor_ca_SET2)
    zkk_SET2_rmsds.append(zkk_rmsd)
    ij8_SET2_rmsds.append(ij8_rmsd)
    apo_SET2_rmsds.append(apo_rmsd)
    inhibitor_SET2_rmsds.append(inhibitor_rmsd)
    
    # Cflank
    zkk_rmsd = md.rmsd(traj_ca_Cflank, zkk_ca_Cflank)
    ij8_rmsd = md.rmsd(traj_ca_Cflank, ij8_ca_Cflank)
    apo_rmsd = md.rmsd(traj_ca_Cflank, apo_ca_Cflank)
    inhibitor_rmsd = md.rmsd(traj_ca_Cflank, inhibitor_ca_Cflank)
    zkk_Cflank_rmsds.append(zkk_rmsd)
    ij8_Cflank_rmsds.append(ij8_rmsd)
    apo_Cflank_rmsds.append(apo_rmsd)
    inhibitor_Cflank_rmsds.append(inhibitor_rmsd)
    
# save results
np.save('CA_rmsds/1zkk_Nflank_rmsds.npy', zkk_Nflank_rmsds)
np.save('CA_rmsds/1zkk_SET1_rmsds.npy', zkk_SET1_rmsds)
np.save('CA_rmsds/1zkk_SETI_rmsds.npy', zkk_SETI_rmsds)
np.save('CA_rmsds/1zkk_SET2_rmsds.npy', zkk_SET2_rmsds)
np.save('CA_rmsds/1zkk_Cflank_rmsds.npy', zkk_Cflank_rmsds)

np.save('CA_rmsds/4ij8_Nflank_rmsds.npy', ij8_Nflank_rmsds)
np.save('CA_rmsds/4ij8_SET1_rmsds.npy', ij8_SET1_rmsds)
np.save('CA_rmsds/4ij8_SETI_rmsds.npy', ij8_SETI_rmsds)
np.save('CA_rmsds/4ij8_SET2_rmsds.npy', ij8_SET2_rmsds)
np.save('CA_rmsds/4ij8_Cflank_rmsds.npy', ij8_Cflank_rmsds)

np.save('CA_rmsds/apo_Nflank_rmsds.npy', apo_Nflank_rmsds)
np.save('CA_rmsds/apo_SET1_rmsds.npy', apo_SET1_rmsds)
np.save('CA_rmsds/apo_SETI_rmsds.npy', apo_SETI_rmsds)
np.save('CA_rmsds/apo_SET2_rmsds.npy', apo_SET2_rmsds)
np.save('CA_rmsds/apo_Cflank_rmsds.npy', apo_Cflank_rmsds)

np.save('CA_rmsds/inhibitor_Nflank_rmsds.npy', inhibitor_Nflank_rmsds)
np.save('CA_rmsds/inhibitor_SET1_rmsds.npy', inhibitor_SET1_rmsds)
np.save('CA_rmsds/inhibitor_SETI_rmsds.npy', inhibitor_SETI_rmsds)
np.save('CA_rmsds/inhibitor_SET2_rmsds.npy', inhibitor_SET2_rmsds)
np.save('CA_rmsds/inhibitor_Cflank_rmsds.npy', inhibitor_Cflank_rmsds)  

# make plots
# 1ZKK
plt.figure()
for i in range(len(zkk_Nflank_rmsds)):
    plt.plot(zkk_Nflank_rmsds[i])
plt.title('CA RMSDs to 1ZKK_Nflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_Nflank_rmsds.pdf')

plt.figure()
for i in range(len(zkk_SET1_rmsds)):
    plt.plot(zkk_SET1_rmsds[i])
plt.title('CA RMSDs to 1ZKK_SET1 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_SET1_rmsds.pdf')

plt.figure()
for i in range(len(zkk_SETI_rmsds)):
    plt.plot(zkk_SETI_rmsds[i])
plt.title('CA RMSDs to 1ZKK_SETI - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_SETI_rmsds.pdf')

plt.figure()
for i in range(len(zkk_SET2_rmsds)):
    plt.plot(zkk_SET2_rmsds[i])
plt.title('CA RMSDs to 1ZKK_SET2 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_SET2_rmsds.pdf')

plt.figure()
for i in range(len(zkk_Cflank_rmsds)):
    plt.plot(zkk_Cflank_rmsds[i])
plt.title('CA RMSDs to 1ZKK_Cflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/1zkk_Cflank_rmsds.pdf')

# 4IJ8
plt.figure()
for i in range(len(ij8_Nflank_rmsds)):
    plt.plot(ij8_Nflank_rmsds[i])
plt.title('CA RMSDs to 4IJ8_Nflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_Nflank_rmsds.pdf')

plt.figure()
for i in range(len(ij8_SET1_rmsds)):
    plt.plot(ij8_SET1_rmsds[i])
plt.title('CA RMSDs to 4IJ8_SET1 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_SET1_rmsds.pdf')

plt.figure()
for i in range(len(ij8_SETI_rmsds)):
    plt.plot(ij8_SETI_rmsds[i])
plt.title('CA RMSDs to 4IJ8_SETI - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_SETI_rmsds.pdf')

plt.figure()
for i in range(len(ij8_SET2_rmsds)):
    plt.plot(ij8_SET2_rmsds[i])
plt.title('CA RMSDs to 4IJ8_SET2 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_SET2_rmsds.pdf')

plt.figure()
for i in range(len(ij8_Cflank_rmsds)):
    plt.plot(ij8_Cflank_rmsds[i])
plt.title('CA RMSDs to 4IJ8_Cflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/4ij8_Cflank_rmsds.pdf')

# Apo
plt.figure()
for i in range(len(apo_Nflank_rmsds)):
    plt.plot(apo_Nflank_rmsds[i])
plt.title('CA RMSDs to Apo_Nflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/apo_Nflank_rmsds.pdf')

plt.figure()
for i in range(len(apo_SET1_rmsds)):
    plt.plot(apo_SET1_rmsds[i])
plt.title('CA RMSDs to Apo_SET1 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/apo_SET1_rmsds.pdf')

plt.figure()
for i in range(len(apo_SETI_rmsds)):
    plt.plot(apo_SETI_rmsds[i])
plt.title('CA RMSDs to Apo_SETI - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/apo_SETI_rmsds.pdf')

plt.figure()
for i in range(len(apo_SET2_rmsds)):
    plt.plot(apo_SET2_rmsds[i])
plt.title('CA RMSDs to Apo_SET2 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/apo_SET2_rmsds.pdf')

plt.figure()
for i in range(len(apo_Cflank_rmsds)):
    plt.plot(apo_Cflank_rmsds[i])
plt.title('CA RMSDs to Apo_Cflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/apo_Cflank_rmsds.pdf')

# Inhibitor
plt.figure()
for i in range(len(inhibitor_Nflank_rmsds)):
    plt.plot(inhibitor_Nflank_rmsds[i])
plt.title('CA RMSDs to Inhibitor_Nflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/inhibitor_Nflank_rmsds.pdf')

plt.figure()
for i in range(len(inhibitor_SET1_rmsds)):
    plt.plot(inhibitor_SET1_rmsds[i])
plt.title('CA RMSDs to Inhibitor_SET1 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/inhibitor_SET1_rmsds.pdf')

plt.figure()
for i in range(len(inhibitor_SETI_rmsds)):
    plt.plot(inhibitor_SETI_rmsds[i])
plt.title('CA RMSDs to Inhibitor_SETI - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/inhibitor_SETI_rmsds.pdf')

plt.figure()
for i in range(len(inhibitor_SET2_rmsds)):
    plt.plot(inhibitor_SET2_rmsds[i])
plt.title('CA RMSDs to Inhibitor_SET2 - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/inhibitor_SET2_rmsds.pdf')

plt.figure()
for i in range(len(inhibitor_Cflank_rmsds)):
    plt.plot(inhibitor_Cflank_rmsds[i])
plt.title('CA RMSDs to Inhibitor_Cflank - 969 trajectories')
plt.xlabel('Frame (0.25 ns / frame)')
plt.ylabel('C-alpha RMSD (nm)')
plt.savefig('CA_rmsds/inhibitor_Cflank_rmsds.pdf')

# now look at the RMSFs
