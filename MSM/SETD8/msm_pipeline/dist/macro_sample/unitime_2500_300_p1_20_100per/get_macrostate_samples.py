# some code copied from choderalab/msm-pipeline

import pyemma
import numpy as np
import mdtraj as md
import glob

project_names = ['SET8_apo']

fnames = glob.glob('../../../../data_cut_start/*/*.h5')

traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_all()

source = pyemma.coordinates.source(fnames, features=feat)

indices = [np.load('samples_300_p1_20_100per.npy')]

for i in range(len(project_names)):
    pyemma.coordinates.save_trajs(source, indices[i], prefix=project_names[i], fmt = 'pdb')
