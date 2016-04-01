import os
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel, mixtape.ghmm
import numpy as np
import mdtraj as md
from parameters import load_trajectories, build_full_featurizer
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils

system = os.path.split(os.getcwd())[1]
out_filename = "/home/kyleb/src/choderalab/pimento/kab/TR01/tex/figures/%s_tics.png" % system

trj0, trajectories, filenames = load_trajectories()

pipeline = sklearn.externals.joblib.load("./pipeline.job")
feature_pipeline = sklearn.externals.joblib.load("./feature_pipeline.job")
cluster_pipeline = sklearn.externals.joblib.load("./cluster_pipeline.job")

X_all = feature_pipeline.transform(trajectories)
q = np.concatenate(X_all)

cluster = cluster_pipeline.steps[-1][1]
covars_ = cluster.covars_
covars_ = cluster.covars_.diagonal(axis1=1, axis2=2)

title("Slow Coordinates in %s" % system)
i, j = 0, 1
hexbin(q[:,i], q[:, j], bins='log')
#errorbar(cluster.means_[:, i], cluster.means_[:, j], xerr=covars_[:,i] ** 0.5, yerr=covars_[:, j] ** 0.5, fmt='kx', linewidth=4)
#offset = np.ones(2) * 0.05
#for state in range(n_states):    
#    plt.annotate("%d" % state, cluster.means_[state, 0:2] + offset, fontsize='x-large')
ylabel("2nd Slowest Coordinate")
xlabel("1st Slowest Coordinate")
savefig(out_filename, bbox_inches="tight")
