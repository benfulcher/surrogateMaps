"""

Generate surrogate expression maps for mouse or human

"""

from src import utils, files, surrogates, spatial
import pandas as pd
import numpy as np

# Parameters:
numMaps = 40000 # generate this many maps
whatSpecies = 'human' # 'mouse', 'human', 'mouseCortex'
whatMode = 'generating' # 'fitting', 'generating'

if whatSpecies=='mouse':
    d0 = 0.78 # length scale (exponential fit from gene coexpression)
    rho = 0.8 # strength of distance effect (fixed)
elif whatSpecies=='mouseCortex':
    d0 = 2.7 # length scale (exponential fit from gene coexpression)
    rho = 0.8 # strength of distance effect (fixed)
elif whatSpecies=='human':
    d0 = 35 # length scale (exponential fit from gene coexpression)
    rho = 0.8 # strength of distance effect (fixed)

# Filenames:
if whatSpecies=='mouse':
    distMatFile = 'mouseDistMat.csv'
    kFile = 'mouseDegree.csv'
    fileNameOut = ('mouseSurrogate_N%u_rho%u_d0%u.csv' % (numMaps,int(rho*10),int(d0*100)))
elif whatSpecies=='mouseCortex':
    distMatFile = 'mouseCortexDistMat.csv'
    kFile = 'mouseCortexDegree.csv'
    fileNameOut = ('mouseCortexSurrogate_N%u_rho%u_d0%u.csv' % (numMaps,int(rho*10),int(d0*100)))
elif whatSpecies=='human':
    distMatFile = 'humanDistMat_99.csv'
    kFile = 'humanDegree_99.csv'
    fileNameOut = ('humanSurrogate_N%u_rho%u_d0%u.csv' % (numMaps,int(rho*10),int(d0*100)))

# Load in the pairwise separation distance matrix:
distMat = pd.read_csv(distMatFile,header=None)
distMat = distMat.to_numpy()
assert distMat.shape[0] == distMat.shape[1]
numAreas = distMat.shape[0]

# Do something:
if whatMode=='fitting':
    # FITTING DEGREE
    k = pd.read_csv(kFile,header=None).values.flatten()
    assert len(k) == numAreas

    rho,d0 = surrogates.fit_parameters(distMat,k)
    print('The loaded phenotype has rho = %f, d0 = %f' % (rho,d0))
else:
    # Generate an ensemble:
    print('Generating %u null maps with d0 = %.2f and rho = %.2f...' % (numMaps,d0,rho))
    nullMaps = np.empty((numAreas,numMaps))
    for i in range(numMaps):
        nullMaps[:,i] = surrogates.generate(distMat,rho,d0)
    # Save out:
    df = pd.DataFrame(nullMaps)
    df.to_csv(fileNameOut)
    print('Saved to %s' % fileNameOut)
