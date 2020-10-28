"""
ePSproc classes selection functions and wrappers

16/10/20    Methods for/from base class.
"""


def Esubset(self, key = None, dataType = None, Erange = None, Etype = None):
    """
    Basic Etype subselection & slice routine for base class and plot wrappers.

    Will return view of data only (http://xarray.pydata.org/en/stable/indexing.html#copies-vs-views).

    TODO: add matEleSelector for thresholding here too?

    Q: currently set for single dataset. May want to handle multiple & set consistent Erange, as per original paradigm.

    Note: originally written for E-subselection, but should work on any slicable dataType + dim.
    NOTE: may have issues using .sel for MultiIndex coord slices if all inds are not specified, see http://xarray.pydata.org/en/stable/indexing.html#multi-level-indexing

    """

    if Erange is None:
        Erange = [self.data[key][dataType][Etype].min().data.item(), self.data[key][dataType][Etype].max().data.item()]

    if len(Erange) == 2:
        Erange.append(None)  # Set step size to None if not passed.

    # More elegant way to swap on dims?
    if Etype == 'Ehv':
    # Subset before plot to avoid errors on empty array selection!
        subset = self.data[key][dataType].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword

    else:
        subset = self.data[key][dataType].sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword

    return subset
