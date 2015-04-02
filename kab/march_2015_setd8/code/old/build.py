import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils
from sklearn.pipeline import make_pipeline
import itertools

STRIDE = 4  # 1 ns
MIN_LENGTH = 400 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

RUN = 0

filenames = [filename for filename in glob.glob(PATH + "run%d-*.h5" % RUN) if len(md.open(filename)) > MIN_LENGTH]
run0 = np.array(map(lambda x: "run1" in x, filenames))
len(filenames)


%time trajectories = [md.load(filename, stride=STRIDE) for filename in filenames]

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")

atom_indices = t_1zkk.top.select("resSeq >=15 and resSeq <= 145 and backbone")
# select ends, not resi 16:144

n_states = 40
n_macrostates = 3

clusterer = cluster.MiniBatchKMedoids(n_clusters=n_states, metric="rmsd")
microstate_model = msm.MarkovStateModel()

pipeline0 = make_pipeline(clusterer, microstate_model)
labels = pipeline0.fit_transform(trajectories)

pcca = lumping.PCCAPlus.from_msm(microstate_model, n_macrostates=n_macrostates)
macrostate_model = msm.MarkovStateModel()
macrostate_labels = macrostate_model.fit_transform(pcca.transform(labels))

sampled_macrostates = macrostate_model.sample_discrete()
selected_pairs_by_state = macrostate_model.draw_samples(macrostate_labels, 5)
sample_trajectories = utils.map_drawn_samples(selected_pairs_by_state, trajectories)

for k, trj in enumerate(sample_trajectories):
    trj.save("state%d.pdb" % k)




import networkx
g = networkx.from_numpy_matrix(macrostate_model.ratemat_, create_using=networkx.DiGraph())

normalize_max = lambda x: x / x.max()
NODE_SCALING = 2000

networkx.draw_spectral(g, node_size=normalize_max(msm_model.populations_) * NODE_SCALING, with_labels=True, font_size=18)
