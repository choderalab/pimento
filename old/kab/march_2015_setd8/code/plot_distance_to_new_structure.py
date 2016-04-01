import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils
from sklearn.pipeline import make_pipeline
import itertools

STRIDE = 4  # 1 ns
MIN_LENGTH = 400 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

RUN = 0
pdb_code = {1:"4ij8", 0:"1zkk"}[RUN]

filenames = [filename for filename in glob.glob(PATH + "run%d-*.h5" % RUN) if len(md.open(filename)) > MIN_LENGTH]
run0 = np.array(map(lambda x: "run1" in x, filenames))
len(filenames)


%time trajectories0 = [md.load(filename, stride=STRIDE) for filename in filenames]

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")
t_new = md.load("../SETD8-DRG_refmac1.pdb")

t_new_bb = t_new.atom_slice(t_new.top.select("protein and backbone and chainid 0"))
t_1zkk_bb = t_1zkk.atom_slice(t_1zkk.top.select("protein and backbone and resSeq >= 2"))
t_4ij8_bb = t_4ij8.atom_slice(t_1zkk.top.select("protein and backbone and resSeq >= 2"))

md.rmsd(t_1zkk_bb, t_new_bb)
md.rmsd(t_4ij8_bb, t_new_bb)
md.rmsd(t_4ij8_bb, t_1zkk_bb)


atom_indices = t_1zkk.top.select("protein and backbone and resSeq >= 2")
# Slide trajectories to align to new file
trajectories = [t.atom_slice(atom_indices) for t in trajectories0]  # For RMSD clustering

r_new = featurizer.RMSDFeaturizer(t_new_bb).transform(trajectories)
r_new_flat = np.concatenate(r_new)
