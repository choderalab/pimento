import mdtraj as md
from msmbuilder import example_datasets, cluster, msm, featurizer, lumping, utils, dataset, decomposition
from sklearn.pipeline import make_pipeline

gamma = 0.09
lag_time = 50

trajectories = dataset.MDTrajDataset("./data/trajectories/run*.h5")
runs = np.array(map(lambda x: "run0" in x, trajectories.glob_matches)).astype('int')
d = dataset.NumpyDirDataset("./data/dihedrals/")

tica = decomposition.tICA(gamma=gamma, n_components=10, lag_time=lag_time)
X = tica.fit_transform(d)

Xf = np.concatenate(X)
Xf0 = np.concatenate(tica.transform([x for k, x in enumerate(d) if runs[k] == 0]))
Xf1 = np.concatenate(tica.transform([x for k, x in enumerate(d) if runs[k] == 1]))

tica.timescales_ 
tica.eigenvalues_

t_1zkk = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/1ZKK_fixed1.pdb")
t_4ij8 = md.load("/home/kyleb/src/choderalab/fah-projects/projects/setd8-p10478/pdbs/4IJ8_fixed1.pdb")
dihedrals_pdb = featurizer.DihedralFeaturizer(types=["phi", "psi", "chi2"]).transform([t_1zkk, t_4ij8])
x_pdb = tica.transform(dihedrals_pdb)


hexbin(Xf[:,0], Xf[:, 1], bins='log')

plot(x_pdb[0][:, 0], x_pdb[0][:, 1], 'o', markersize=20, label="1zkk", color='purple')
plot(x_pdb[1][:, 0], x_pdb[1][:, 1], 'o', markersize=20, label="4ij8", color='white')

title("Dihedral tICA Analysis")
xlabel("Slowest Coordinate")
ylabel("Second Slowest Coordinate")

legend(loc=0)

savefig("./figures/dihedral_map.png", bbox_inches="tight")


figure()
hexbin(Xf[:,0], Xf[:, 1], bins='log')

plot(Xf0[:, 0], Xf0[:, 1], 'x', color="black", label="run 0")
plot(Xf1[:, 0], Xf1[:, 1], 'x', color="gray", label="run 1")

plot(x_pdb[0][:, 0], x_pdb[0][:, 1], 'o', markersize=20, label="1zkk", color='purple')
plot(x_pdb[1][:, 0], x_pdb[1][:, 1], 'o', markersize=20, label="4ij8", color='white')

title("Dihedral tICA Analysis")
xlabel("Slowest Coordinate")
ylabel("Second Slowest Coordinate")

legend(loc=0)

savefig("./figures/dihedral_map_byrun.png", bbox_inches="tight")
