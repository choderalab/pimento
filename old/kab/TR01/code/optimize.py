import glob
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel
import numpy as np
import mdtraj as md
from mixtape import ghmm, feature_selection, subset_featurizer
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils

n_iter = 1000
n_rounds = 100
n_choose = 50
stride = 1
lag_time = 1
n_components = 2

filenames = glob.glob("./Trajectories/*.h5")
trajectories = [md.load(filename) for filename in filenames]

if len(trajectories) > 1:
    trajectories = trajectories[0::2]
else:
    t = trajectories[0]
    trajectories = [trajectories[0][0:trajectories[0].n_frames/2]]

try:
    featurizer = sklearn.externals.joblib.load("./featurizer-%d-%d.job" % (n_components, n_choose))
except:
    featurizer = mixtape.subset_featurizer.guess_featurizers(trajectories[0][0], n_choose)

model = mixtape.tica.tICA(lag_time=lag_time, n_components=n_components)
tica_optimizer = mixtape.feature_selection.Optimizer(featurizer, model, n_iter)

for i in range(n_rounds):
    print("i = %d" % i)
    featurizer = tica_optimizer.optimize(trajectories)
    print("saving %d" % i)
    sklearn.externals.joblib.dump(featurizer, "./featurizer-%d-%d.job" % (n_components, n_choose), compress=True)
