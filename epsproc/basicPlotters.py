"""
ePSproc basic plotting functions

Some basic functions for 2D/3D plots.

07/11/19    v1


TODO
----


"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Basic plotting for molecular structure
# TODO: replace with more sophisticated methods!
def molPlot(molInfo):
    """Basic 3D scatter plot from molInfo data."""

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4], s = 100* molInfo['atomTable'][:,0])
    ax.plot(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4])
    ax.axis('off');

    # Add labels
    lShift = 0.01* np.array([1., 1., 1.])
    for row in molInfo['atomTable']:
        ax.text(row[2] + lShift[0], row[3] + lShift[1], row[4] + lShift[2], f"Z = {row[0]}") # molInfo['atomTable'][:,0])

    plt.show()
