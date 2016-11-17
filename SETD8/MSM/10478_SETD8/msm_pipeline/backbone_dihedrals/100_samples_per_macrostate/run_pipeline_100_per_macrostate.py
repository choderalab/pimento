from msmpipeline import run_pipeline
import glob

fnames = glob.glob('../../data/*.h5')

run_pipeline(fnames, project_name='setd8', n_structures_per_macrostate=100)
