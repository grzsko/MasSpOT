# coding: utf-8

# In[4]:


with open("data/ptm/seq.txt") as infile:
    fasta = ""
    for line in infile:
        fasta += line.strip()


# In[5]:


from masserstein.peptides import get_protein_counter
from masserstein import Spectrum

from collections import Counter
import random
import itertools
from copy import copy

import gc

from optimal_transport_py import quick_distance


# In[6]:


htt_counter = get_protein_counter(fasta)
htt_spectrum = Spectrum.new_from_fasta(fasta, threshold=0.1)
htt_mass = htt_spectrum.average_mass()


# In[7]:


changes = {
    "S": Counter({"H": -1, "P": 1, "O": 3}),
    "K": Counter({"H": 2, "C": 2, "O": 1}),
    "T": Counter({"H": -1, "P": 1, "O": 3})
}

changes_counts = {
    "S": 30,
    "K": 11
}

ids = {
    0: "S",
    1: "K"
}


# In[8]:


count = 100



randoms = random.sample(list(itertools.product(*(range(i + 1) for i in changes_counts.values()))), count)


# In[9]:


mass_dists = []
spec_dists = []

gc.collect()

def do_work(arg):
    i, t = arg
    new_counter = copy(htt_counter)
    for index, coor in enumerate(t):
        change = changes[ids[index]]
#         print(change, new_counter)
        for _ in range(coor):
            new_counter += change
    new_formula = ''.join(sym+str(count) for sym,count in sorted(new_counter.items()))
    new_spectrum = Spectrum(new_formula, threshold=0.1)

    print("Doing:", i)
    spec_dist = quick_distance(htt_spectrum, new_spectrum)
    mass_dist = new_spectrum.average_mass() - htt_mass
    gc.collect()
    return spec_dist, mass_dist


if __name__ == "__main__":
    # from multiprocessing import Pool
    # pool = Pool(20)
    for sd, md in map(do_work, enumerate(randoms)):
        mass_dists.append(md)
        spec_dists.append(sd)


# In[10]:


    import json
    with open("data/ptm/results.json", "w") as outfile:
        json.dump(outfile, [spec_dists, mass_dists])


# In[8]:




# In[ ]:




