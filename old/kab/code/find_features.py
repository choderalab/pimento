import mixtape.featurizer, mixtape.tica
import numpy as np
import mdtraj as md
from mixtape import ghmm, selector, subset_featurizer, selector
from parameters import load_trajectories, build_full_featurizer
import sklearn.pipeline, sklearn.externals.joblib

n_iter = 5000
n_choose = 40
lag_time = 4

trj0, trajectories, filenames = load_trajectories()

featurizer = build_full_featurizer(trj0, n_choose)
featurizer = sklearn.externals.joblib.load("./featurizer_full.job")

tica_optimizer = mixtape.selector.TICAOptimizer(featurizer, trajectories, lag_time=lag_time)
tica_optimizer.optimize(n_iter, trajectories)

featurizer = tica_optimizer.featurizer
sklearn.externals.joblib.dump(featurizer, "./featurizer_full.job", compress=True)
