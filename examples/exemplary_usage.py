#!/usr/bin/env python
# coding: utf-8

# In[1]:


# get_ipython().run_line_magic('load_ext', 'autoreload')
#
# get_ipython().run_line_magic('autoreload', '2')


# In[2]:


# from pyimzml.ImzMLParser import ImzMLParser
from masserstein import Spectrum

# from sklearn.cluster import AgglomerativeClustering
# from sklearn.metrics import pairwise_distances
import os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

from entropic_regularized_distance import calculate_distance

from datetime import datetime

from optimal_transport_py import quick_distance

# import matplotlib.pyplot as plt


# In[3]:


# Get data

# p = ImzMLParser('data/mouse_cebellum/test_POS.imzML')
# dimensions = (max(coor[0] for coor in p.coordinates),
#               max(coor[1] for coor in p.coordinates))
# # We know that z has only one value - 1
# picture = [[None for y in range(dimensions[1])] for x in range(dimensions[0])]
# for idx, (x, y, z) in enumerate(p.coordinates):
#     mzs, intensities = p.getspectrum(idx)
#     s = Spectrum('', empty=True, label=str(x) + ", " + str(y))
#     s.confs = list(zip(mzs, intensities))
#     picture[x - 1][y - 1] = s
dimensions = (81, 21)
picture = [[None for y in range(dimensions[1])] for x in range(dimensions[0])]
with open("data/mouse_cebellum/peak_picked_brain.tsv", "r") as infile:
    mzs = list(map(float, next(infile).strip().split("\t")[2:]))
    for line in infile:
        line = line.strip().split("\t")
        x, y = map(int, line[:2])
        intensities = list(map(float, line[2:]))
        s = Spectrum('', empty=True, label=str(x) + ", " + str(y))
        s.confs = [x for x in zip(mzs, intensities) if x[1] > 0]
        s.normalize()
        picture[x - 1][y - 1] = s


# In[42]:


# dimensions


# In[4]:


def flatten(l):
    return [item for sublist in l for item in sublist]


# In[5]:


flat_picture = flatten(picture)

# def own_affinity():
#     count = 0
#     distances = np.zeros((len(flat_picture), len(flat_picture)))
#     for i, row in enumerate(distances):
#         for j in range(i + 1):
#             gc.collect()
#             distances[i][j] = calculate_distance(flat_picture[i],
#                                                  flat_picture[j])
#             count += 1
#             if count % 100 == 0:
#                 print(datetime.datetime.now(), count, i, j)
#     return distances



# distance_matrix = own_affinity()



# In[54]:


# gc.collect()
# # from entropic_regularized_distance import global_vs
# import cProfile
# cProfile.run('calculate_distance(flat_picture[80], flat_picture[20])')
# # calculate_distance(flat_picture[4], flat_picture[3])


# In[8]:





def par_distance(arg):
    i, j = arg
    res = []
    return i, j, quick_distance(flat_picture[i], flat_picture[j])



def ind_generator():
    count = 0
    for i in range(len(flat_picture)):
        for j in range(len(flat_picture)):
            count += 1
            if count % 1000 == 0:
                print("Now yielding ", count)
            yield i, j
    return


def own_affinity_paralel():
    import gc
    gc.collect()
    from multiprocessing import Pool
    pool = Pool(64)
    distance_matrix = \
        [[None for y in range(len(flat_picture))] for x in range(len(flat_picture))]
    for i, j, res in pool.imap_unordered(par_distance, ind_generator(),
                                         chunksize=100):
        distance_matrix[i][j] = res
    return distance_matrix


def own_affinity():
    import gc
    gc.collect()
    distance_matrix = \
        [[None for y in range(len(flat_picture))] for x in range(len(flat_picture))]
    for i, j, res in map(par_distance, ind_generator()):
        distance_matrix[i][j] = res
    return distance_matrix


if __name__ == "__main__":
    start = datetime.now()
    distance_matrix = own_affinity()

    print("Calculated, took", datetime.now() - start)
    # In[9]:


    # import json
    # with open("data/mouse_cebellum/distances.json", "w") as outfile:
    #     json.dump(distance_matrix, outfile)


# In[15]:


# # print(res)
# # flat_picture[4].plot()
# # flat_picture[3].plot()
# get_ipython().run_line_magic('matplotlib', 'notebook')
# # plt.plot(range(len(global_vs)), global_vs)
# flat_picture[80].plot()
#
#
# # In[75]:
#
#
# count = 0
# distances = np.zeros((len(flat_picture), len(flat_picture)))
# for i, row in enumerate(distances):
#         for j in range(i + 1):
#             count += 1
#             print(i, j,
#                   len(flat_picture[i].confs), len(flat_picture[j].confs),
#                   count)
#             if count == 13:
#                 raise Exception()
#
#
#
# # In[ ]:
#
#
# ac = AgglomerativeClustering(n_clusters=None, affinity='precomputed',
#                              linkage='average')
# ac.fit(distance_matrix)
#
#
# # In[ ]:
#
#
#
#
#
# # In[19]:
#
#
# import gc
# gc.collect()
#
#
# # In[ ]:
#
#
#
#
