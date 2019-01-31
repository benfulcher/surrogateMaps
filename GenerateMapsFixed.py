"""

Generate surrogate expression maps for mouse or human

"""

from src import utils, files, surrogates, spatial
import pandas as pd
import numpy as np

# Parameters:
whatSpecies = 'human' # 'mouse', 'human'
if whatSpecies=='mouse':
    d0 = 0.4 # length scale
    rho = 0.8 # strength of distance effect
elif whatSpecies=='human':
    d0 = 20 # length scale
    rho = 0.8 # strength of distance effect

numMaps = 10000 # generate this many maps

# File names:
if whatSpecies=='mouse':
    distMatFile = 'mouseDistMat.csv'
    fileNameOut = ('mouseSurrogate_rho%u_d0%u.csv' % (int(rho*10),int(d0*100)))
elif whatSpecies=='human':
    distMatFile = 'humanDistMat.csv'
    fileNameOut = ('humanSurrogate_rho%u_d0%u.csv' % (int(rho*10),int(d0*100)))

# Load in the pairwise separation distance matrix:
distMat = pd.read_csv(distMatFile,header=None)
assert distMat.shape[0] == distMat.shape[1]
numAreas = distMat.shape[0]

# Generate
print('Generating %u null maps with d0 = %.2f and rho = %.2f...' % (numMaps,d0,rho))
nullMaps = np.empty((numAreas,numMaps))
for i in range(numMaps):
    nullMaps[:,i] = surrogates.generate(distMat,rho,d0)

# Save out:
df = pd.DataFrame(nullMaps)
df.to_csv(fileNameOut)
print('Saved to %s' % fileNameOut)
