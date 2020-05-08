from MasSpOT import spectral_distance
import numpy as np
from multiprocessing import Pool
import gc
from sklearn.cluster import AgglomerativeClustering

def own_affinity(flat_picture):
    count = 0
    distances = np.zeros((len(flat_picture), len(flat_picture)))
    for i, row in enumerate(distances):
        for j in range(i + 1):
            gc.collect()
            distances[i][j] = spectral_distance(flat_picture[i], flat_picture[j])
            count += 1
            # if count % 100 == 0:
            # print(datetime.datetime.now(), count, i, j)
    return distances


def flatten(l):
    return [item for sublist in l for item in sublist]


# def par_distance(arg):
#     i, j = arg
#     return spectral_distance(flat_picture[i], flat_picture[j]), i, j


def ind_generator(flat_picture):
    count = 0
    for i in range(len(flat_picture)):
        for j in range(len(flat_picture)):
            count += 1
            if count % 10 == 0:
                print("Now yielding ", count)
                return
            yield (i, j)
    return

#
# def own_affinity_paralel():
#     pool = Pool(2)
#     distance_matrix = [[None for y in range(len(flat_picture))] for x in range(len(flat_picture))]
#     for dist, i, j in pool.imap_unordered(par_distance, ind_generator(), chunksize=2):
#         distance_matrix[i][j] = dist
#     return distance_matrix


#%%
def perform_clusterization(picture, dimensions):
    flat_picture = flatten(picture)

    distance_matrix = own_affinity(flat_picture)
    ac = AgglomerativeClustering(n_clusters=None, affinity='precomputed',
                                 linkage='average', compute_full_tree=True,
                                 distance_threshold=0.45)
    labels = ac.fit_predict(distance_matrix)

    label_picture = [[None for y in range(dimensions[1])] for x in range(dimensions[0])]

    indexes = [(i, j) for i, sublist in enumerate(picture) for j, item in enumerate(sublist)]

    for n, (i, j) in zip(labels, indexes):
        label_picture[i][j] = n

    return label_picture