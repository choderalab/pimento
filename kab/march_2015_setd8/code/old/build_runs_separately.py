import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils
from sklearn.pipeline import make_pipeline
import itertools

STRIDE = 4  # 1 ns
MIN_LENGTH = 400 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

RUN = 1
pdb_code = {1:"4ij8", 0:"1zkk"}[RUN]

filenames = [filename for filename in glob.glob(PATH + "run%d-*.h5" % RUN) if len(md.open(filename)) > MIN_LENGTH]
run0 = np.array(map(lambda x: "run1" in x, filenames))
len(filenames)


%time trajectories0 = [md.load(filename, stride=STRIDE) for filename in filenames]

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")

atom_indices = t_1zkk.top.select("resSeq >=15 and resSeq <= 145")
# select ends, not resi 16:144
trajectories = [t.atom_slice(atom_indices) for t in trajectories0]  # For RMSD clustering


#atraj = trajectories0[0][0].join([t[::150] for t in trajectories0])
#atraj.save("./all_run%d.pdb" % RUN)

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
sample_trajectories = utils.map_drawn_samples(selected_pairs_by_state, trajectories0)

for k, trj in enumerate(sample_trajectories):
    trj.save("run_%d-state%d.pdb" % (RUN, k))
