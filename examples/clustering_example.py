
from masserstein import Spectrum

from pyimzml.ImzMLParser import ImzMLParser
#%%

# # Obtain data
# !mkdir msi_data
# !wget -O msi_data/test_POS.imzML https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS487/download/9c756a3f-2c96-4449-8dd7-d64540df5c6c\?file\=test_POS.imzML
# !wget wget -O msi_data/test_POS.ibd https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS487/download/9c756a3f-2c96-4449-8dd7-d64540df5c6c\?file\=test_POS.ibd

#%%
# Parse data
p = ImzMLParser('msi_data/test_POS.imzML')
dimensions = (max(coor[0] for coor in p.coordinates),
              max(coor[1] for coor in p.coordinates))
# We know that z has only one value: 1
picture = [[None for y in range(dimensions[1])] for x in range(dimensions[0])]
for idx, (x, y, z) in enumerate(p.coordinates):
    mzs, intensities = p.getspectrum(idx)
    s = Spectrum(confs=list(zip(mzs, intensities)), label=str(x - 1) + ", " + str(y - 1))
    # remove peptide artifacts
    s.confs = [x for x in s.confs if x[0] < 1000]
    picture[x - 1][y - 1] = s

#%%
# Apply peak-picking procedure
for row in picture:
    for spectrum in row:
        spectrum.confs = spectrum.find_peaks()
        spectrum.confs = spectrum.centroid(0.5)

#%%
from MasSpOT import perform_clusterization
label_picture = perform_clusterization(picture, dimensions)

#%%
from matplotlib import pyplot as plt
fig, ax = plt.subplots()
# x =
im = ax.imshow(label_picture, cmap=plt.get_cmap("Oranges"))

fig.colorbar(im)

# fig.tight_layout()
plt.show()