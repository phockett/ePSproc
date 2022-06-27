"""
ePSproc additional plotting utility functions

13/01/22

"""

# Show plot from tmo-dev - assumes dict based data structure
# def showPlot(self, dim, run, pType='curve'):
#     """Render (show) specified plot from pre-set plot dictionaries (very basic!)."""
#
#     try:
#         if self.__notebook__:
#             display(self.data[run][pType][dim])  # If notebook, use display to push plot.
#         else:
#             return self.data[run][pType][dim]  # Otherwise return hv object.
#
#     except:
#         print(f'Plot [{run}][{pType}][{dim}] not set.')

# Basic show plot routine
def showPlot(plotObj, returnPlot = False, __notebook__ = False):
    """Render (show) specified plot if notebook (very basic!)."""

    # try:
    if __notebook__:
        display(plotObj)  # If notebook, use display to push plot.

    if returnPlot or not __notebook__:
        return plotObj  # Return hv object.

    # except:
    #     print(f'Plot [{run}][{pType}][{dim}] not set.')
