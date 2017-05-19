import pyemma
import msmbuilder
from msmbuilder import decomposition
import msmpipeline
from msmpipeline.contact_features import find_respairs_that_changed
import numpy as np
import glob
import mdtraj as md

fnames = glob.glob('../data_cut_start/*/*.h5')
traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)

scheme = 'closest'
respairs_that_changed = find_respairs_that_changed(fnames, scheme=scheme)
np.save('respairs_that_changed', respairs_that_changed)

feat.add_residue_mindist(residue_pairs=respairs_that_changed, scheme=scheme)
source = pyemma.coordinates.source(fnames, features = feat)

X = source.get_output()

# sparse tICA - default rho
sparse_tica = msmbuilder.decomposition.SparseTICA(n_components=50, lag_time=50, kinetic_mapping=True)
sparse_tica.fit(X)

sparse_eigenv = sparse_tica.eigenvectors_
np.save('sparse_eigenv', sparse_eigenv)
sparse_eigenvl = sparse_tica.eigenvalues_
np.save('sparse_eigenvl', sparse_eigenvl)

sparse_Y = sparse_tica.transform(X)
np.save('sparse_Y', sparse_Y)
del sparse_Y

# sparse tICA - two orders of magnitude lower rho (1e-4)
sparse_tica_lowrho = msmbuilder.decomposition.SparseTICA(n_components=50, lag_time=50, kinetic_mapping=True, rho=0.0001)
sparse_tica_lowrho.fit(X)

sparse_lowrho_eigenv = sparse_tica_lowrho.eigenvectors_
np.save('sparse_lowrho_eigenv', sparse_lowrho_eigenv)
sparse_lowrho_eigenvl = sparse_tica_lowrho.eigenvalues_
np.save('sparse_lowrho_eigenvl', sparse_lowrho_eigenvl)

sparse_lowrho_Y = sparse_tica_lowrho.transform(X)
np.save('sparse_lowrho_Y', sparse_lowrho_Y)
del sparse_lowrho_Y

# normal tICA
tica = msmbuilder.decomposition.tICA(n_components=50, lag_time=50, kinetic_mapping=True)
tica.fit(X)

eigenv = tica.eigenvectors_
np.save('normal_tica_eigenv', eigenv)
eigenvl = tica.eigenvalues_
np.save('normal_tica_eigenvl', eigenvl)

Y = tica.transform(X)
np.save('normal_tica_Y', Y)
del Y

# get starting conf coordinates
# with no hydrogen models
fnames = ['starting_structures/all_traj.h5']
traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_residue_mindist(residue_pairs=respairs_that_changed, scheme=scheme)

start_conf_source = pyemma.coordinates.source(fnames, features=feat)
start_conf_X = start_conf_source.get_output()

start_conf_Y_sparse_tica = sparse_tica.transform(start_conf_X)
start_conf_Y_sparse_tica_lowrho = sparse_tica_lowrho.transform(start_conf_X)
start_conf_Y_normal_tica = tica.transform(start_conf_X)

np.save('start_conf_Y_sparse_tica_noH_models', start_conf_Y_sparse_tica)
np.save('start_conf_Y_sparse_tica_lowrho_noH_models', start_conf_Y_sparse_tica_lowrho)
np.save('start_conf_Y_normal_tica_noH_models', start_conf_Y_normal_tica)

# with hydrogen containing models - openmm.app.Modeller hydrogenated
fnames = ['starting_structures/all_traj_H_modeller.h5']
traj = md.load(fnames[0])
top = traj.top

feat = pyemma.coordinates.featurizer(top)
feat.add_residue_mindist(residue_pairs=respairs_that_changed, scheme=scheme)

start_conf_source = pyemma.coordinates.source(fnames, features=feat)
start_conf_X = start_conf_source.get_output()

start_conf_Y_sparse_tica = sparse_tica.transform(start_conf_X)
start_conf_Y_sparse_tica_lowrho = sparse_tica_lowrho.transform(start_conf_X)
start_conf_Y_normal_tica = tica.transform(start_conf_X)

np.save('start_conf_Y_sparse_tica_H_modeller', start_conf_Y_sparse_tica)
np.save('start_conf_Y_sparse_tica_lowrho_H_modeller', start_conf_Y_sparse_tica_lowrho)
np.save('start_conf_Y_normal_tica_H_modeller', start_conf_Y_normal_tica)

# training scores
scores = [sparse_tica.score_, sparse_tica_lowrho.score_, tica.score_]
np.save('gmrq_scores', scores)

# save X - in portions to not overblow memory
np.save('X[:500]', X[:500])
np.save('X[500:1000]', X[500:1000])
np.save('X[1000:1500]', X[1000:1500])
np.save('X[1500:2000]', X[1500:2000])
np.save('X[2000:2500]', X[2000:2500])                                                                                            
np.save('X[2500:3000]', X[2500:3000])
np.save('X[3000:3500]', X[3000:3500])
np.save('X[3500:4000]', X[3500:4000])
np.save('X[4000:4500]', X[4000:4500])
np.save('X[4500:]', X[4500:])
