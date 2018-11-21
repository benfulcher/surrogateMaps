"""

Generating surrogate maps for mouse expression patterns

"""

from src import utils, files, surrogates, spatial
import pandas as pd
import numpy as np


# Parameters:
d0_range = (0.01,0.1,0.2,0.4,0.6,0.8,1,10,100,1000) # length scales
rho = 0.5 # effect size
numMaps = len(d0_range)

# File names:
distMatFile = 'mouseDistMat.csv'
fileNameOut = ('mouseSurrogate_rho%u.csv' % int(rho*10))


# Load in the pairwise separation distance matrix:
distMat = pd.read_csv(distMatFile,header=None)
assert distMat.shape[0] == distMat.shape[1]
numAreas = distMat.shape[0]

# Generate
print("# Generating %u surrogate maps across a range of rho values... " % numMaps)
nullMaps = np.empty((numMaps, numAreas))
for i,d0 in enumerate(d0_range):
    print('d0 = %g' % d0)
    nullMaps[i,:] = surrogates.generate(distMat, rho, d0)

# Save out:
df = pd.DataFrame(nullMaps,index=d0_range)
df.to_csv(fileNameOut)
print('Saved to %s' % fileNameOut)
