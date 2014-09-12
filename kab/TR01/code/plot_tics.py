import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel, mixtape.ghmm
import numpy as np
import mdtraj as md
from parameters import load_trajectories, build_full_featurizer
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils

system = "NSD2"
out_filename = "/home/kyleb/src/choderalab/pimento/kab/TR01/tex/figures/%s_tics.png" % system

trj0, trajectories, filenames = load_trajectories()


X_all = feature_pipeline.transform(trajectories)
q = np.concatenate(X_all)

covars_ = cluster.covars_
covars_ = cluster.covars_.diagonal(axis1=1, axis2=2)

for i, j in [(0, 1)]:
    figure()
    title("%d" % n_states)
    hexbin(q[:,i], q[:, j], bins='log')
    errorbar(cluster.means_[:, i], cluster.means_[:, j], xerr=covars_[:,i] ** 0.5, yerr=covars_[:, j] ** 0.5, fmt='kx', linewidth=4)

offset = np.ones(2) * 0.05
for state in range(n_states):    
    plt.annotate("%d" % state, cluster.means_[state, 0:2] + offset, fontsize='x-large')



pipeline = sklearn.externals.joblib.load("./partial_pipeline.job")
full_pipeline = sklearn.externals.joblib.load("./full_pipeline.job")

X_all = pipeline.fit_transform(trajectories)
q = np.concatenate(X_all)

model = full_pipeline.steps[-1][1]

title("Slow Coordinates in %s" % system)
i, j = 0, 1
hexbin(q[:,i], q[:, j], bins='log')
ylabel("2nd Slowest Coordinate")
xlabel("1st Slowest Coordinate")
savefig(out_filename, bbox_inches="tight")




#figure()
#i, j = 0, 1
#hexbin(q[:,i], q[:, j], bins='log')

#offset = np.ones(2) * 0.05
#errorbar(model.means_[states, i], model.means_[states, j], xerr=model.vars_[states,i] ** 0.5, yerr=model.vars_[states, j] ** 0.5, fmt='kx', linewidth=4)


"""
i, j = 0, 1
hexbin(q[:,i], q[:, j], bins='log')
errorbar(model.means_[:, i], model.means_[:, j], xerr=model.vars_[:,i] ** 0.5, yerr=model.vars_[:, j] ** 0.5, fmt='kx', linewidth=4)

offset = np.ones(2) * 0.05
for state in range(model.n_states):
    plt.annotate("%d" % state, model.means_[state, 0:2] + offset, fontsize='x-large')
"""
