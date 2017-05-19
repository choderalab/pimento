import pyemma
import numpy as np

Y = list(np.load('../../setd8_distances_tica_projection.npy'))

clustering = pyemma.coordinates.cluster_uniform_time(data=Y, k=2500)

np.save('clustercenters', clustering.clustercenters)

dtrajs = [dtraj.flatten() for dtraj in clustering.get_output()]

np.save('dtrajs', dtrajs)


