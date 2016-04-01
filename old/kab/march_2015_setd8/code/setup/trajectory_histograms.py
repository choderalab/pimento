import glob
import mdtraj as md

PATH = "/home/kyleb/dat/fah_data/10478/"

RUN = 0

x = np.array([len(md.open(filename)) for filename in glob.glob(PATH + "run*.h5")]) / 4.

sum(x[x>300])
mean(x[x>300])
sum(x>300)
