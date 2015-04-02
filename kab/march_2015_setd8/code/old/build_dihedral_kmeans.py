import pandas as pd
import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils, decomposition
from sklearn.pipeline import make_pipeline
import itertools

STRIDE = 4  # 1 ns
MIN_LENGTH = 400 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

RUN = 1
pdb_code = {1:"4ij8", 0:"1zkk"}[RUN]

filenames = [filename for filename in glob.glob(PATH + "run%d-*.h5" % RUN) if len(md.open(filename)) > MIN_LENGTH]
%time trajectories = [md.load(filename, stride=STRIDE) for filename in filenames]

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")



n_macrostates = 3


dfeaturizer = featurizer.DihedralFeaturizer()
dihedrals = dfeaturizer.transform(trajectories)

gamma = 1000.
lag_time = 1
tica = decomposition.tICA(n_components=5, lag_time=lag_time, gamma=gamma)
X = tica.fit_transform(dihedrals)
Xflat = np.concatenate(X)
hexbin(Xflat[:, 0], Xflat[:, 1])
tica.timescales_

n_states = 10
clusterer = cluster.KMeans(n_clusters=n_states)
microstate_model = msm.MarkovStateModel(lag_time=lag_time)
p2 = make_pipeline(clusterer, microstate_model)
sequences = p2.fit_transform(X)
microstate_model.score(sequences)


f = lambda x, k: sum(x[0:k + 1]) / sum(x)

gridvals = range(0, 10)

for timestep in [1, 5, 10, 50, 100, 200]:
    gmrq_ratio = pd.Series(np.array([f(microstate_model.eigenvalues_ ** timestep, k) for k in gridvals]), index=gridvals)
    gmrq_ratio.plot(style='o', label=timestep)

legend()
