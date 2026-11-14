# This script demonstrates uses of MTEX for EBSD data corrections.
# Run inside a virtual environment with the required packages installed.
# Requires MTEX, H5py and a GUI backend of your choosing to run. Other imports are optional.
# Requires python 3.12.x - 3.14.x. This program has been tested using Python 3.13.x

# Script is based of the MATLAB version written by Thomas Bennett and Monica Martinez Henriquez
# Adapted to Python and revised by Joseph Quarshie, for use in the Taleff Group 2024.06.12

# Imports
import glob # Finding data files
import matplotlib
import sys
matplotlib.use = lambda *a, **k: None # Prevent MTEX from changing the backend
matplotlib.use('QtAgg') # use QtAgg or macosx. Use Agg for headless rendering.
from matplotlib.backends.qt_compat import QtWidgets
import matplotlib.pyplot as plt
import numpy as np
from GrainMerge import * # Custom function for selecting grains.
from mtex import * # Lazy import. Can fix for performance.
from pytictoc import TicToc # Measuring execution time

# Prevent last window closing from exiting program
if QtWidgets.QApplication.instance() is None:
    QtWidgets.QApplication.instance() == QtWidgets.QApplication(sys.argv)
QtWidgets.QApplication.instance().setQuitOnLastWindowClosed(False)

t = TicToc() #create instance of class

## Using relative paths, specify the location of the EBSD data file relative to script.
dir = 'EBSD_Files/Sample'
dataFiles = glob.glob(f'{dir}/**/*.ang', recursive=True)
matName = 'Molybdenum'

# Define the crystal structure and symmetry
CS = crystalFrame('m-3m', [2.87, 2.87, 2.87], mineral = matName, color='LightSkyBlue')

# Define the specimen symmetry
SS = specimenFrame('orthorhombic') # pyMTEX seemes to have dropped this system

# Set the orientation convention for plotting
# ARROWS_UTF8 = ("←", "→", "↑", "↓", "⊗", "⊙")
# ARROWS_ASCII = ("<-", "->", "^", "v", "(x)", "(.)")
plottingConvention.default('y^->x') # Adjust based on EBSD data orientation

## Load in EBSD Data
# Import EBSD data and make the appropriate correction
ebsdR = EBSD.load(dataFiles[0], CS, interface='ang')
print(ebsdR)

# Convert EBSD data to a grid to make it easier to work with
ebsdR = gridify(ebsdR)

# Define the inverse pole figure key
ipfkey = ipfHSVKey(ebsdR[matName])
ipfkey.ipfDirection = zvector
rgb = ipfkey.orientation2color(ebsdR[matName].orientations)

# Plot the IPF key
ipfkey.plot()
plt.gcf().canvas.manager.set_window_title(f'IPF Key - {matName}')

# Save IPF Key for reference
plt.savefig(f'{dir}/{matName}_IPF_key.png')

# Generate the IPF map
plot(ebsdR[matName], rgb)
plt.gcf().canvas.manager.set_window_title(f'IPF Map - {matName}')
# Save IPF Map for reference
plt.savefig(f'{dir}/{matName}_IPFMap_raw.png')
plt.close()

# IQ Plot Normalization for better visualization
iq = ebsdR[matName].iq
iqsort = np.sort(iq)
iqThres = iqsort[int(np.ceil(.99 * len(iqsort)))]
iq[np.isnan(iq)] = 0
iq = (iq - np.min(iq)) / (iqThres - np.min(iq))
iq[iq > 1] = 1

# Image Quality Plot
im = plot(ebsdR[matName], num2rgb(iq, colormap = black2whiteColorMap), cmap='gray')
plt.gcf().colorbar(im, ax=plt.gca(), label='IQ')   
plt.gcf().canvas.manager.set_window_title(f'IQ Map - {matName}')
# Save corrected IQ map for reference
plt.savefig(f'{dir}/{matName}_IQ_corrected.png')

# Confidence Index Plot, setting to black and white
ci = ebsdR[matName].ci
im = plot(ebsdR[matName], num2rgb(ci, colormap = black2whiteColorMap), cmap='gray')
plt.gcf().colorbar(im, ax=plt.gca(), label='CI')
plt.gcf().canvas.manager.set_window_title(f'CI Map - {matName}')
# Save corrected CI map for reference
plt.savefig(f'{dir}/{matName}_CI_corrected.png')

# Field Signal Detector Plot, setting to black and white
im = plot(ebsdR, num2rgb(ebsdR.sem_signal, colormap = black2whiteColorMap), cmap='gray')
plt.gcf().colorbar(im, ax=plt.gca(), label='FSD')
plt.gcf().canvas.manager.set_window_title(f'FSD Map - {matName}')
# Save FSD map for reference
plt.savefig(f'{dir}/{matName}_FSD.png')

# Generate IPF Key and IQ Plot 
rgb_ipf = rgb * iq[..., np.newaxis]  # darken by IQ
plot(ebsdR[matName], rgb_ipf)
plt.gcf().canvas.manager.set_window_title(f'IPF + IQ - {matName}')
plt.savefig(f'{dir}/{matName}_IPF_IQ.png')

# Thresholding data by confidence index
CI = ebsdR.ci
minCI = 0.1

# Preserve ebsdR before tampering 
ebsd = ebsdR.copy()
ebsd.phaseId[CI < minCI] = 1 # Alter low-confidence points to the single phase.

# Set requirements to count as a grain from the raw data, removing small grains
grains = calcGrains(ebsd, angle=5 * degree, minPixel=10) # investigate this data type

# Main loop for correcting mis-indexed grains. Turn interactive mode back on
#plt.ion()
while True:
    plt.close('all')
    plot(ebsd[matName], ebsd[matName].orientations, micronbar='on')
    hold(True)
    largeGrains = grains[grains.numPixel > 20]
    text(largeGrains, largeGrains.id)
    plot(grains.boundary, lineWidth=1)
    children, parents = select_grains(grains, matName) # Pull children and parent grains out to combine them.
    print(children, parents) # We will just call the function again if user wants to make more selections.
    pairs = grains.neighbors() # This is important for validating user selections.

    idx = find_boundary_index(pairs, children, parents) # Finds specific boundary between selected grains.
    print(f"\nFound at index {idx}\n" if idx != -1 else "\n Not found\n")

    # Find the largest selected parent grain and grab its mean orientation
    L_parent = 0
    L_parent_meanO = None
    for i in range(len(parents)):
        parent_area = float(np.asarray(grains['id', parents[i]].area).item())
        if parent_area > L_parent:
            L_parent = parent_area
            L_parent_meanO = grains['id', parents[i]].meanOrientation

    ind = np.isin(ebsd.grainId, children) # Index child grains
    ebsd.rotations.data[ind] = L_parent_meanO.data # Set children orientations to the parent's mean orientation.
    grains = calcGrains(ebsd, angle=5 * degree, minPixel=10) # Recalculate grains

    if input("Continue? [y/n] ").lower() != 'y':
        plt.close('all')
        break


# Smooth grains (altering data) using a half quadratic filter
t.tic()
quad = halfQuadraticFilter()
quad.alpha = 0.5
ebsd = smooth(ebsd, quad, fill=True)
t.toc('\nApplying a half quadratic filter took')

# Fill missing points once more to correct points near boundaries
t.tic()
ebsd = fill(ebsd, extrapolate=True)
t.toc('Filling missing points took')

# Determine and smooth grain boundaries again
t.tic()
grains = calcGrains(ebsd, angle=5 * degree, minPixel=10)
grains = smoothBoundary(grains,10)
t.toc('Determining and smoothing grain boundaries took')

# Generate the IPF map
colors = ipfkey.orientation2color(ebsd.orientations)
plot(ebsd, colors, micronbar='off')
hold(True)
plot(grains.boundary, lineWidth=1)
hold(False)
plt.gcf().canvas.manager.set_window_title(f'Corrected IPF Map - {matName}')
# Save the corrected IPF map
plt.savefig(f'{dir}/{matName}_IPFMap_Corrected.png')

plot(ebsd[matName], ipfkey.orientation2color(ebsd[matName].orientations), micronbar='on')
plt.gcf().canvas.manager.set_window_title(f'Smoothed Grains - {matName}')
hold(True)
plot(grains.boundary, lineWidth=2)
hold(False)

# Plot the grain reference orientation deviation (GROD) map
t.tic()
grod = calcGROD(ebsd, grains)
t.toc('Calculating GROD took')
im = plot(ebsd, angle(grod) / degree)
plt.gcf().colorbar(im, ax=plt.gca(), label='GROD')
hold(True)
plot(grains.boundary, lineWidth=1)
hold(False)
plt.gcf().canvas.manager.set_window_title(f'GROD Map - {matName}')
plt.savefig(f'{dir}/{matName}_GROD.png')

# Plot the kernel average misorientation (KAM) map
t.tic()
kam = ebsd.KAM() / degree
t.toc('Calculating KAM took')

print(f'\nDefault KAM median: {np.nanmedian(kam):.2f} degree; maximum: {np.nanmax(kam):.1f} degree\n')

plot(ebsd, kam, micronbar=False)
plt.gcf().canvas.manager.set_window_title(f'KAM Map - {matName}')
setColorRange([0, 15])
mtexColorMap('LaboTeX')
mtexColorbar(title='KAM in degree')
hold(True)
plot(grains.boundary, lineWidth=1.5)
hold(False)
# Save corrected KAM map for reference
plt.savefig(f'{dir}/{matName}_KAM_corrected.png')

# Plot final grain boundaries
n = len(grains.boundary) 
plot(grains.boundary, lineWidth=1.5, color=[(0.541, 0.169, 0.882)] * n)    # must have a color for each boundary segment
plt.gcf().canvas.manager.set_window_title(f'Final Grain Boundaries - {matName}')
plt.savefig(f'{dir}/{matName}_Final_Grain_Boundaries.png')

# Set a function to estimate the orientation distribution function (ODF)
t.tic()
psi = calcKernel(grains.meanOrientation)
t.toc('Determining the optimal kernel halfwidth took')
print(f'\nThe optimal kernel halfwidth is {psi.halfwidth()/degree} degrees.\n')
OptAng = psi.halfwidth() / degree

# Apply the optimal kernel halfwidth to the ODF estimation
t.tic()
odf = calcDensity(ebsd.orientations, halfwidth=OptAng * degree)
t.toc('Calculating ODF took')

# Determine Miller indices for pole figure generation
h = Miller([1,0,0], [1,1,0], [1,1,1], odf.CS)   
plotPF(odf, h, antipodal=True, smooth=True, upper=True)
#mtexColorMap('parula')
mtexColorbar()
plt.gcf().canvas.manager.set_window_title(f'Pole Figures - {matName}')
plt.savefig(f'{dir}/{matName}_Pole_Figures.png')

print('All done! All maps have been generated and saved.')