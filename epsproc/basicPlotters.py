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

# For Arrow3D class
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# Arrow3D class from https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
# Code: https://stackoverflow.com/a/22867877
# Builds on existing 2D arrow class, projects into 3D.
# From StackOverflow user CT Zhu, https://stackoverflow.com/users/2487184/ct-zhu
class Arrow3D(FancyArrowPatch):
    """
    Define Arrow3D plotting class

    Code verbatim from StackOverflow post https://stackoverflow.com/a/22867877
    Thanks to CT Zhu, https://stackoverflow.com/users/2487184/ct-zhu

    """

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# Basic plotting for molecular structure
# TODO: replace with more sophisticated methods!
def molPlot(molInfo):
    """Basic 3D scatter plot from molInfo data."""

    #*** Plot atoms
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4], s = 100* molInfo['atomTable'][:,0], c = molInfo['atomTable'][:,0])
    ax.axis('off');

    # Add labels + bonds for small systems (need more logic for larger systems!)
    if molInfo['atomTable'].shape[0] < 4:
        ax.plot(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4])

        lShift = 0.01* np.array([1., 1., 1.])
        for row in molInfo['atomTable']:
            ax.text(row[2] + lShift[0], row[3] + lShift[1], row[4] + lShift[2], f"Z = {row[0]}") # molInfo['atomTable'][:,0])

    else:
        # Generate legend from scatter data, as per https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/scatter_with_legend.html
        legend1 = ax.legend(*scatter.legend_elements(),loc="lower left", title="Z")
        ax.add_artist(legend1)


    #*** Plot axes with direction vecs

    dirs = ('x','y','z')  # Labels
    oV = np.zeros([3,3])   # Generate origin
    xyzV = np.identity(3)  # Generate unit vectors

    # Scale for lines, check axis limits & molecular extent
    # sf = 1
    axLims = np.asarray([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
    sfAxis = np.asarray([ax.get_xlim()[1], ax.get_ylim()[1], ax.get_zlim()[1]])  # Get +ve axis limits
    sfMol = 1.8* np.asarray([molInfo['atomTable'][:,2].max(), molInfo['atomTable'][:,3].max(), molInfo['atomTable'][:,4].max()]) # Get max atom position

    # Set sensible limits
    mask = sfMol > 0
    sfPlot = sfMol*mask + (sfAxis/2)*(~mask)

    for n, xyz in enumerate(xyzV):
        sf = sfPlot[n]
    #     ax.plot([oV[n,0], sf*xyz[0]], [oV[n,1], sf*xyz[1]], [oV[n,2], sf*xyz[2]])  # Basic line
    #     ax.quiver(oV[n,0],oV[n,1],oV[n,2],sf*xyz[0],sf*xyz[1],sf*xyz[2], arrow_length_ratio = 0.05, lw=2)  # Use quiver for arrows - not very pleasing as arrow heads suck
        a = Arrow3D([oV[n,0], sf*xyz[0]], [oV[n,1], sf*xyz[1]], [oV[n,2], sf*xyz[2]], mutation_scale=10,
                    lw=4, arrowstyle="-|>", color="r", alpha=0.5)  # With Arrow3D, defined above
        ax.add_artist(a)
        ax.text(sf*xyz[0], sf*xyz[1], sf*xyz[2], dirs[n])


    # Add (y,z) plane as a surface plot
    yy, zz = np.meshgrid(range(-1, 2), range(-1, 2))
    xx = np.zeros([3,3])
    ax.plot_surface(xx*sfPlot[0],yy*sfPlot[1],zz*sfPlot[2], alpha=0.1)

    # Reset plot limits (otherwise rescales to arrows!)
    ax.set_xlim(axLims[0,:])
    ax.set_ylim(axLims[1,:])
    ax.set_zlim(axLims[2,:])

    # TO CONSIDER:
    # See https://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D
    # ax.get_proj()  # Get projection matrix
    # ax.view_init(0,20)  # Set view, (az,el)

    plt.show()
