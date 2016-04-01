import glob
import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils
from sklearn.pipeline import make_pipeline
import itertools

STRIDE = 4  # 1 ns
MIN_LENGTH = 400 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

filenames = [filename for filename in glob.glob(PATH + "*.h5") if len(md.open(filename)) > MIN_LENGTH]
run0 = np.array(map(lambda x: "run1" in x, filenames))
len(filenames)


%time trajectories = [md.load(filename, stride=STRIDE) for filename in filenames]

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")

atom_indices = t_1zkk.top.select("resSeq >=15 and resSeq <= 145 and backbone")
# select ends, not resi 16:144


r_1zkk = featurizer.RMSDFeaturizer(t_1zkk, atom_indices=atom_indices).transform(trajectories)
r_4ij8 = featurizer.RMSDFeaturizer(t_4ij8, atom_indices=atom_indices).transform(trajectories)

r_1zkk_flat = np.concatenate(r_1zkk)
r_4ij8_flat = np.concatenate(r_4ij8)


for VALUE in [True, False]:
    figure()
    pdb_code = dict(True="4ij8", False="1zkk")[str(VALUE)]
    hexbin(r_1zkk_flat, r_4ij8_flat)
    for k in range(len(filenames)):
        x = r_1zkk[k]
        y = r_4ij8[k]
        if run0[k] == VALUE:
            plot(x[:, 0], y[:, 0], 'kx')
    xlabel("RMSD to 1zkk [nm]")
    ylabel("RMSD to 4ij8 [nm]")
    title("RMSD Distribution starting from %s" % pdb_code)
    savefig("./figures/rmsd_from_%s.png" % pdb_code, bbox_inches="tight")
