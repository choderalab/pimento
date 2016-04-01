import glob
import mdtraj as md
import os

MIN_LENGTH = 300 * 4
PATH = "/home/kyleb/dat/fah_data/10478/"

for RUN in [0, 1]:
    filenames = [filename for filename in glob.glob(PATH + "run%d-*.h5" % RUN) if len(md.open(filename)) > MIN_LENGTH]
    for filename in filenames:
        base_filename = os.path.split(filename)[1]
        out_filename = "./data/run%d/%s" % (RUN, base_filename)
        if not os.path.exists(out_filename):
            os.symlink(filename, out_filename)


filenames = [filename for filename in glob.glob(PATH + "run*.h5") if len(md.open(filename)) > MIN_LENGTH]
for filename in filenames:
    base_filename = os.path.split(filename)[1]
    out_filename = "./data/trajectories/%s" % (base_filename)
    if not os.path.exists(out_filename):
        os.symlink(filename, out_filename)
