import glob
import mixtape.featurizer, mixtape.tica, mixtape.cluster, mixtape.markovstatemodel, mixtape.ghmm
import numpy as np
import mdtraj as md
import sklearn.pipeline, sklearn.externals.joblib
import mixtape.utils

n_choose = 50
stride = 1
lag_time = 1
n_components = 2

filenames = glob.glob("./Trajectories/*.h5")
trajectories = [md.load(filename) for filename in filenames]

featurizer = sklearn.externals.joblib.load("./featurizer-%d-%d.job" % (n_components, n_choose))

n_states = 7
tica = mixtape.tica.tICA(n_components=n_components, lag_time=lag_time)
subsampler = mixtape.utils.Subsampler(lag_time=lag_time)
msm = mixtape.markovstatemodel.MarkovStateModel(n_timescales=5)
cluster = mixtape.cluster.GMM(n_components=n_states, covariance_type='full')
feature_pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica)])
cluster_pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica), ("cluster", cluster)])
pipeline = sklearn.pipeline.Pipeline([("features", featurizer), ('tica', tica), ("subsampler", subsampler), ("cluster", cluster), ("msm", msm)])

pipeline.fit(trajectories)

sklearn.externals.joblib.dump(pipeline, "./pipeline.job", compress=True)
sklearn.externals.joblib.dump(cluster_pipeline, "./cluster_pipeline.job", compress=True)
sklearn.externals.joblib.dump(feature_pipeline, "./feature_pipeline.job", compress=True)


states = cluster_pipeline.transform(trajectories)
ind = msm.draw_samples(states, 3)
samples = mixtape.utils.map_drawn_samples(ind, trajectories)

for i in range(n_states):
    for k, t in enumerate(samples[i]):
        t.save("pdbs/state%d-%d.pdb" % (i, k))
