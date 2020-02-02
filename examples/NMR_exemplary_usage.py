#!/usr/bin/env python
# coding: utf-8

# In[34]:


# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')
# 
# import matplotlib.pyplot as plt

# get_ipython().run_line_magic('matplotlib', 'notebook')

import statistics


# In[18]:


import os
import numpy as np
from pygam import LinearGAM, s

from masserstein import Spectrum


# In[3]:


# Load spectra

def load_spectra(prefix):
    result = []
    for filename in os.listdir(prefix):
        with open(prefix + filename, "r") as infile:
            confs = []
            next(infile) # header
            for line in infile:
                line = line.strip()
                if line != "":
                    line = line.split()
                    line = list(map(float, line))
                    # zakresOD, zakresDO, calka
                    if line[1] != 0:
                        confs.append((line[0], line[1]))
        s = Spectrum(None, empty=True, label=filename)
        s.confs = confs
        result.append(s)
    return result
        
drugs = load_spectra("data/nmr/data/ASCI/leki/")
active = load_spectra("data/nmr/data/ASCI/czynne/")
nonactive = load_spectra("data/nmr/data/ASCI/pomocnicze/")


# In[4]:


# one_confs = drugs[0].confs
# diffs = [(v[0], v[1] - one_confs[i-1][1]) for i, v in enumerate(one_confs[1:])]
# diffed_spec = Spectrum(formula=None, empty=True)
# diffed_spec.confs = diffs


# In[5]:


# #diffed_spec.plot(profile=True)

# abs_diffs = np.abs(list(zip(*diffs))[1])
# ppm_s, signal = list(zip(*one_confs))
# not_signal_idx = np.where(abs_diffs <= np.quantile(np.abs(list(zip(*diffs))[1]), 0.9))[0]
# not_signal_idx = np.where(signal <= np.quantile(signal, 0.85))[0]
# print(not_signal_idx)

# noise_spec = Spectrum(formula=None, empty=True)
# noise_spec.confs = [one_confs[i] for i in not_signal_idx]
# noise_spec.plot(profile=True)
# ppm_noise, signal_noise = zip(*noise_spec.confs)
# # print(ppm_noise, signal_noise)
# # Spectrum.plot_all([noise_spec, drugs[0]], profile=True)


# In[6]:




# gam = LinearGAM(s(0))
# gam.gridsearch(np.array([[i] for i in ppm_noise]), np.asarray(signal_noise))


# In[7]:


# gam.summary()
# prediction = gam.predict(X=np.array([[i] for i in ppm_s]))


# In[8]:


# baseline_spectrum = Spectrum(formula=None, empty=True, label="cokolwiek")
# baseline_spectrum.confs = [(ppm_s[i], -v) for i, v in enumerate(prediction)]
# # result_spectrum.plot(profile=True)

# result_spectrum = drugs[0] + baseline_spectrum
# result_spectrum.plot(profile=True)
# # Spectrum.plot_all([baseline_spectrum, drugs[0]], profile=True)


# In[9]:


# result_spectrum.normalize()
# result_spectrum.cut_smallest_peaks(0.1)


# In[10]:


# result_spectrum.plot(profile=True)


# In[11]:


def find_baseline(spectrum, save=False, show=False):
    one_confs = spectrum.confs
    diffs = [(v[0], v[1] - one_confs[i - 1][1])
             for i, v in enumerate(one_confs[1:])]
    diffed_spec = Spectrum(formula=None, empty=True)
    diffed_spec.confs = diffs
    # abs_diffs = np.abs(list(zip(*diffs))[1])
    ppm_s, signal = list(zip(*one_confs))
    not_signal_idx = np.where(signal <= np.quantile(signal, 0.85))[0]

    noise_spec = Spectrum(formula=None, empty=True)
    noise_spec.confs = [one_confs[i] for i in not_signal_idx]
    ppm_noise, signal_noise = zip(*noise_spec.confs)
    gam = LinearGAM(s(0))
    gam.gridsearch(np.array([[i] for i in ppm_noise]),
                   np.asarray(signal_noise))
    prediction = gam.predict(X=np.array([[i] for i in ppm_s]))
    baseline_spectrum = Spectrum(formula=None, empty=True, label="Baseline")
    baseline_spectrum.confs = [(ppm_s[i], v) for i, v in enumerate(prediction)]

    result_spectrum = spectrum + -1 * baseline_spectrum

    Spectrum.plot_all([baseline_spectrum, spectrum], profile=True, show=show)
    if save:
        plt.savefig("/vagrant/pictures/pic" + spectrum.label + ".pdf")
    return baseline_spectrum, result_spectrum
    
def denoise(spectrum):
    spectrum.normalize()
    spectrum.cut_smallest_peaks(0.1)
    return spectrum


# In[106]:


all_spectra = drugs + active + nonactive


def remove_pure_noises(spectra):
    res = []
    inds = []
    
    for i, sp in enumerate(spectra):
        value = abs(max(c[1] for c in sp.confs) / min(c[1] for c in sp.confs))
        if not value < 2:
            res.append(sp)
            inds.append(i)
    return res, inds

# all_spectra, all_inds = remove_pure_noises(all_spectra)


# In[107]:


def find_baselines(spectra):
    baselines = []
    debaselined = []
    for sp in spectra:
        b, r = find_baseline(sp, show=False, save=False)
        baselines.append(b)
        debaselined.append(r)
    return baselines, debaselined

def find_baselines_p(spectra):
    baselines = []
    debaselined = []
    from multiprocessing import Pool
    pool = Pool(3)
    for b, r in pool.imap(find_baseline, spectra, chunksize=10):
        # b, r = find_baseline(sp, show=False, save=False)
        baselines.append(b)
        debaselined.append(r)
    return baselines, debaselined


all_used_spectra = []
import json
with open("data/nmr/data/debaselined.json", "r") as infile:
    data = json.load(infile)
    for d in data:
        sp = Spectrum(formula=None, empty=True, label=d[0])
        sp.confs = d[1]
        all_used_spectra.append(sp)
# all_baselines, all_spectra = find_baselines(all_spectra)
for sp in all_used_spectra:
    sp.normalize()

# In[112]:


# one_s= all_spectra[125]
# one_s.plot()

# from copy import deepcopy
# second_s = deepcopy(one_s)
# print(len(second_s.confs))
# second_s.normalize()
# second_s.cut_smallest_peaks(0.001)
# print(len(second_s.confs))
# second_s.plot()

# for i, sp in enumerate(drugs):
#     value= sum(c[1] for c in sp.confs) / len(sp.confs)
#     if value < 20:
#         print(i, value)
# iter = nonactive
# for i, sp in enumerate(iter):
# #     if sp.label == "Mannitol2.txt":
# #         sp.plot()
# #         print(sum(c[1] for c in sp.confs) / len(sp.confs))
# #         print(max(c[1] for c in sp.confs), min(c[1] for c in sp.confs))
# #         break
#     value = abs(max(c[1] for c in sp.confs) / min(c[1] for c in sp.confs))
#     if value < 2:
#         print(i, value)

# iter[449].plot()
# sum([1 for sp in all_spectra if sum([1 for c in sp.confs if c[1] < 0]) > 0])


# In[123]:


from optimal_transport_py import quick_distance

def par_distance(arg):
    i, j = arg
    res = []
    return i, j, quick_distance(all_used_spectra[i], all_used_spectra[j])


def ind_generator():
    count = 0
    for i in range(len(all_used_spectra)):
        for j in range(len(all_used_spectra)):
            count += 1
            if count % 1000 == 0:
                print("Now yielding ", count)
            yield i, j
    return


def own_affinity_paralel():
    import gc
    gc.collect()
    from multiprocessing import Pool
    pool = Pool(70)
    distance_matrix = [[None for y in range(len(all_used_spectra))] for x in range(len(all_used_spectra))]
    for i, j, res in pool.imap_unordered(par_distance, ind_generator(),
                                         chunksize=1000):
        distance_matrix[i][j] = res
    return distance_matrix


# In[124]:



if __name__ == "__main__":
    from datetime import datetime
    start = datetime.now()
    distance_matrix = own_affinity_paralel()
    print("Calculated, took", datetime.now() - start)
    import json
    with open("data/nmr/data/distances_normalized.json", "w") as outfile:
        json.dump(distance_matrix, outfile)


# In[114]:


# import json
# get_ipython().run_line_magic('pinfo', 'json.dump')


# In[ ]:




