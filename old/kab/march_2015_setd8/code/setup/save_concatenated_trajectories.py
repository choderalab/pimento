import glob
import mdtraj as md

stride = 750
stride2 = 4

for run in [0, 1]:
    filenames = glob.glob("./data/run%d/*.h5" % run)
    traj = md.load(filenames, stride=stride)
    traj[::stride2].save("./run%d.pdb" % run)
