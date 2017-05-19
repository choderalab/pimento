from msmpipeline import run_pipeline
import glob

fnames = glob.glob('../../data_cut_start/*/*.h5')

run_pipeline(fnames, project_name='setd8_distances', n_structures_per_macrostate=100, feature_selection='residue-mindist')
