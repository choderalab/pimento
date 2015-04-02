import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils, dataset, decomposition
from sklearn.pipeline import make_pipeline

runflag = "bothruns"
#gamma = 0.09
gamma = 0.0325 # run1
lag_time = 50

d = dataset.NumpyDirDataset("./data/dihedrals-%s-withchi2/" % runflag)

tica = decomposition.tICA(gamma=gamma, n_components=10, lag_time=lag_time)
X = tica.fit_transform(d)


Xf = np.concatenate(X)

tica.timescales_ 
tica.eigenvalues_


clusterer = cluster.GMM(n_components=3)
sequences = clusterer.fit_transform(map(lambda x: x[:, 0:2], X))

m = msm.MarkovStateModel(lag_time=lag_time)
m.fit(sequences)


