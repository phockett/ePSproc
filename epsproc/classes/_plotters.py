"""
ePSproc classes plot function wrappers

16/10/20    Methods for/from base class.

TODO: should be able to simplify with subselection logic removed (bit ad hoc currently due to organic growth of this!) and/or rewrite as decorators.


"""

from matplotlib import pyplot as plt  # For plot legends with Xarray plotter
import numpy as np # Needed only for np.nan at the moment

from epsproc import matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX, sphFromBLMPlot
from epsproc import lmPlot as lmPlotCore  # Hack rename here to prevent circular logic with local function - TODO: fix with core fn. reorg.
from epsproc.util.env import isnotebook

# ************** Plotters

def plotGetCro(self, pType = 'SIGMA', Erange = None, Etype = 'Eke', keys = None, backend = 'mpl'):
    """
    Basic GetCro (cross-section) data plotting for multijob class. Run self.plot.line(x=Etype, col='Type') for each dataset.
    (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCroComp()` for comparitive plots over datasets.)

    Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.

    Parameters
    ----------
    pType : str, optional, default = 'SIGMA'
        Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).
        If backend = 'hv' this parameter is not used.

    Erange : list of int or float, optional, default = None
        Set plot range [Emin, Emax]. Defaults to full data range if not set

    Etype : str, optional, default = 'Eke'
        Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

    keys : list, optional, default = None
        Keys for datasets to plot.
        If None, all datasets will be plotted.

    backend : str, optional, default = 'mpl'
        Set plotter to use.

        - 'mpl' : Use Matplotlib/native Xarray plotter
        - 'hv' : use Holoviews via :py:func:`epsproc.plotters.hvPlotters.XCplot()`

    """
    # Default to all datasets
    if keys is None:
        keys = self.data.keys()

#         if self.jobs['jobStructure'] == 'subDirs':
    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
        # for m, item in enumerate(self.data[key]['XS']):

        # Set default to full range, same for all cases
        if Erange is None:
            Erange = [self.data[key]['XS'][Etype].min().data, self.data[key]['XS'][Etype].max().data]

        subset = self.data[key]['XS']
        jobLabel = self.data[key]['XS'].jobLabel

        # More elegant way to swap on dims?
        if Etype == 'Ehv':
        # Subset before plot to avoid errors on empty array selection!
            subset = subset.swap_dims({'Eke':'Ehv'}) #.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

            # if subset.any():
            #     subset.plot.line(x=Etype, col='Type')

        subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        if subset.any():
            # Propagate attrs for plot labelling
            subset.attrs = self.data[key]['XS'].attrs

            if backend == 'mpl':
                subset.sel(XC=pType).plot.line(x=Etype, col='Type')

                if pType == 'SIGMA':
                    plt.ylabel('XS/Mb')
                else:
                    plt.ylabel(r"$\beta_{LM}$")

            if backend == 'hv':
                layout, *_ = hvPlotters.XCplot(subset, kdims = Etype)

                print(f"\n*** {jobLabel}")
                # Check if notebook and output
                if isnotebook():
                    display(layout) # Use IPython display(), works nicely for notebook output
                    # Q: is display() always loaded automatically? I think so.
                else:
                    layout  # Not yet tested - not sure what the options are here for hv/Bokeh.
                            # May just want to return this object (hv layout)?

                # Setting title here doesn't work, now set in XCplot()
                # display(layout.opts(title=jobLabel))

#         else:




# Basically as per plotGetCro, but subselect and put on single plot.
# Should do all of this with Holoviews...!
def plotGetCroComp(self, pType='SIGMA', pGauge='L', pSym=('All','All'), Erange = None, Etype = 'Eke',  keys = None):
    """
    Basic GetCro (cross-section) data plotting for multijob class, comparitive plots.
    Run self.plot.line(x=Etype) for each dataset after subselection on Gauge and Symmetry, and use single axis.
    (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCro()` for plots per dataset.)

    Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.

    Parameters
    ----------
    pType : str, optional, default = 'SIGMA'
        Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).

    pGauge : str, optional, default = 'L'
        Set gauge, either 'L' (Length), 'V' (Velocity) or 'M' (Mixed)

    pSym : tuple of strs, optional, default = ('All','All')
        Select symmetry, (Cont, Targ).
        Default value will plot all allowed symmetries.

    Erange : list of int or float, optional, default = None
        Set plot range [Emin, Emax]. Defaults to full data range if not set

    Etype : str, optional, default = 'Eke'
        Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

    keys : list, optional, default = None
        Keys for datasets to plot.
        If None, all datasets will be plotted.

    """



    # Comparison plots over orbs

    # from matplotlib import pyplot as plt  # For legend
    lText = []

    # Default to all datasets
    if keys is None:
        keys = self.data.keys()

    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
        # for m, item in enumerate(self.data[key]['XS']):

        # Set default to full range, same for all cases
        if Erange is None:
            Erange = [self.data[key]['XS'][Etype].min().data, self.data[key]['XS'][Etype].max().data]


        # More elegant way to swap on dims?
        if Etype == 'Ehv':
            # Subset before plot to avoid errors on empty array selection!
            subset = self.data[key]['XS'].swap_dims({'Eke':'Ehv'}).sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword

            # if subset.any():
            #     pltObj = subset.plot.line(x=Etype)

        else:
            subset = self.data[key]['XS'].sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        if subset.any():
            pltObj = subset.plot.line(x=Etype)
            lText.append(self.data[key]['jobNotes']['orbLabel'])

        # Label with orb_sym
#                 lText.append(self.dataSets[key]['XS'][m].attrs['fileBase'].rsplit('/',maxsplit=1)[0])

        # lText.append(f"Orb {self.dataSets[key]['XS'][m].attrs['orbInfo']['orbN']} ({self.dataSets[key]['XS'][m].attrs['orbInfo']['orbSym'][0]})")

        # lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])

    plt.legend(lText)

    if pType == 'SIGMA':
        plt.ylabel('XS/Mb')
    else:
        plt.ylabel(r"$\beta_{LM}$")




def lmPlot(self, Erange = None, Etype = 'Eke', dataType = 'matE', xDim = None, keys = None, refDataKey = None, reindexTol = 0.5, reindexFill = np.nan, setPD = True, **kwargs):
    """
    Wrapper for :py:func:`epsproc.lmPlot` for multijob class. Runs lmPlot() for each dataset.

    Parameters
    ----------
    Erange : list of int or float, optional, default = None
        Set plot range [Emin, Emax]. Defaults to full data range if not set

    Etype : str, optional, default = 'Eke'
        Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

    dataType : str, optional, default = 'matE'
        Set data type to plot, corresponding to label in self.data
        - 'matE' raw matrix elements.
        - 'AFBLM' computed AF BLMs.

    xDim : str, optional, default = None
        Settings for x-axis, if None plot vs. Etype.
        See :py:func:`epsproc.lmPlot` for more details.

    keys : list, optional, default = None
        Keys for datasets to plot.
        If None, all datasets will be plotted.

    refDataKey : str or int, optional, default = None
        If set, calculate difference plots against reference dataset.
        This must be a key in self.data.
        TODO: implement difference plots.
        TODO: implement testing logic, may fail without E-axis forcing, and sym summation?

    reindexTol : float, optional, default = 0.1
        If computing difference data, the reference data is reindexed to ensure E grid matching.
        This specifies tolerance (in E units, usually eV) for reindexing.
        If this fails, difference plot may be null.

    reindexFill : int or float, optional, default = NaN
        Value to use for missing values upon reindexing.
        Default matches [Xarray.reindex default](http://xarray.pydata.org/en/stable/generated/xarray.DataArray.reindex.html), i.e. NaN, but this may give issues in some cases.

    setPD : bool, optional, default = True
        Set Pandas array in main dataset?

    kwargs : dict, optional, default = {}
        Plotting options to pass to :py:func:`epsproc.lmPlot`.
        These will also be set in self.lmPlotOpts for further use.
        Note that any existing options in self.lmPlotOpts will also be used, or overwritten if matching keys are found.


    Notes
    -----
    Basic scheme from ePSmultijob.plotGetCro, which loops and switches on Eke/Ehv. Should tidy up at some point.

    """
    # Set xDim if not passed.
    if xDim is None:
        xDim = Etype

    # Default to all datasets
    if keys is None:
        keys = list(self.data.keys())

    # Set lmPlotOpts
    # Check passed args vs. self.lmPlotOpts and overwrite
    if kwargs:
        for key, value in kwargs.items():
            self.lmPlotOpts[key] = value

    # Set default to full range of 1st dataset, keep same for all cases
    # TODO: set per dataset?
    if Erange is None:
        Erange = [self.data[keys[0]][dataType][Etype].min().data, self.data[keys[0]][dataType][Etype].max().data]

    # Set ref dataset if required
    if refDataKey is not None:
        refData = self.data[refDataKey][dataType]

        if Etype == 'Ehv':
            refData = refData.swap_dims({'Eke':'Ehv'})

        refData = refData.sel(**{Etype:slice(Erange[0], Erange[1])})
        refData.attrs = self.data[refDataKey][dataType].attrs # Propagate atrrs.

    else:
        refData = None


    # Loop over datasets
    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works

        # Init empty list for daPlotpd data
        if setPD:
            self.data[key]['daPlotpd'] = []

        # for m, item in enumerate(self.data[key]['matE']):

        # More elegant way to swap on dims?
        if Etype == 'Ehv':
        # Subset before plot to avoid errors on empty array selection!
            subset = self.data[key][dataType].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        else:
            subset = self.data[key][dataType].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        # Difference data
        # NOTE: ref. data is reindexed here to ensure E-point subtraction (otherwise mismatch will give null results)
        # TODO: may want some selection logic here, this may produce empty results in many cases.
        if refData is not None:
            subset = subset - refData.reindex(**{Etype:subset[Etype]}, method='nearest', tolerance=reindexTol, fill_value=reindexFill)
            subset.attrs = self.data[key][dataType].attrs  # Propagate attribs
            subset.attrs['jobLabel'] = f"{subset.attrs['jobLabel']} diff with {refData.attrs['jobLabel']}"

            # Slightly crude check for empty result.
            # Note this is different from subset.any() below, which only checks for 0
            # See https://numpy.org/doc/stable/reference/generated/numpy.any.html
            if subset.max().isnull():
                print("*** Warning, difference array is Null. Try a larger reindexTol value.")

        # Run lmPlot
        if subset.any():
            daPlot, daPlotpd, legendList, gFig = lmPlotCore(subset, xDim = xDim, **self.lmPlotOpts)
        else:
            daPlotpd = None

        # Set Pandas table to dataset if specified.
        if setPD:
            self.data[key]['daPlotpd'].append(daPlotpd)  # Set to include None cases to keep indexing. Should set as dict instead?



def BLMplot(self, Erange = None, Etype = 'Eke', dataType = 'AFBLM', xDim = None, keys = None, ):
    """
    Basic BLM line plots using Xarray plotter.

    See https://epsproc.readthedocs.io/en/latest/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html

    Similar to :py:func:`epsproc.BLMplot`, but without thresholding, and no dim selection.

    For more flexibility, use `self.lmPlot`.

    TODO: update BLMplot to support more datatypes, and implement here instead.

    TODO: fix dim handling and subselection, see old plotting code.

    """
    # Set xDim if not passed.
    if xDim is None:
        xDim = Etype

    # Default to all datasets
    if keys is None:
        keys = list(self.data.keys())

    # Set default to full range of 1st dataset, keep same for all cases
    # TODO: set per dataset? Now in subset fn.
    # if Erange is None:
    #     Erange = [self.data[keys[0]][dataType][Etype].min().data, self.data[keys[0]][dataType][Etype].max().data]

    # Loop over datasets
    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works

        # for m, item in enumerate(self.data[key]['matE']):

        # # More elegant way to swap on dims?
        # if Etype == 'Ehv':
        # # Subset before plot to avoid errors on empty array selection!
        #     subset = self.data[key][dataType].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
        #
        # else:
        #     subset = self.data[key][dataType].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        subset = self.Esubset(key = key, dataType = dataType, Etype = Etype, Erange = Erange)

        if subset.any():
            if hasattr(subset, 'XSrescaled'):
                print(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, XS")
                subset.XSrescaled.real.plot(x=Etype, col='Labels')
                # plt.title(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, XS")

            print(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, {dataType}")
            subset.real.plot(x=Etype, col='Labels', row='BLM')  # Nice... should give line plots or surfaces depending on dims
            # plt.title(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, {dataType}")


# Plot PADs from mat elements (MF) or BLMs
def padPlot(self, selDims = {}, sumDims = {}, Erange = None, Etype = 'Eke',
                keys = None, dataType = 'TX',
                pType = 'a2', pStyle = 'polar',
                backend = 'mpl'):

    """
    Plot PADs from mat elements (MF) or BLMs.
    """

    # Default to all datasets
    if keys is None:
        keys = self.data.keys()
    else:
        if not isinstance(keys, list):   # Force list if single item passed
            keys = [keys]

    # Loop over datasets
    for key in keys:
        # for m, item in enumerate(self.dataSets[key]['TX']):
        plotFlag = True

        # # More elegant way to swap on dims?
        # subset = self.data[key][dataType]
        #
        # if Etype == 'Ehv':
        # # Subset before plot to avoid errors on empty array selection!
        #     subset = subset.swap_dims({'Eke':'Ehv'})
        #
        # # Slice on Erange with/without step size - bit ugly.
        # if Erange is not None:
        #     if len(Erange) == 3:
        #         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword
        #     else:
        #         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})
        subset = self.Esubset(key = key, dataType = dataType,  Erange = Erange, Etype = Etype)

        # Handling for Euler or Labels dim
        eDim = 'Euler'
        if hasattr(subset, 'Labels'):
#                 if 'Labels' in selDims:
            if eDim in subset.dims:
                subset = subset.swap_dims({'Euler':'Labels'})

            eDim = 'Labels'

        # if selDims is not None:  # Now set as empty dict to avoid issues later.
        # TODO: smarter dim handling here - selDims may not be consistent over all dataSets!
        if selDims:
            subset = subset.sel(selDims)

        # Check dimensionality - sum over Sym & it if necessary
        if (subset.ndim > 3):  # Need some additional dim logic here!
            print(f"Found dims {subset.dims}, summing to reduce for plot. Pass selDims to avoid.")
            if 'Sym' in subset.dims:
                subset = subset.sum('Sym').squeeze()
            if 'it' in subset.dims:
                subset = subset.sum('it').squeeze()

            # Check result is OK
            if (subset.ndim > 3) and not (eDim in subset.dims):
                print(f"*** ERROR: dataset {self.data[key][dataType].jobLabel}  ndims = {subset.ndim}, dims = {subset.dims}. Skipping MFPAD plotting.")
                plotFlag = False

            # else:
            #     pass

        # Propagate attrs
        subset.attrs = self.data[key][dataType].attrs

        # Compute PADs if not already present in dataset
        if dataType != "TX":
            subset, _ = sphFromBLMPlot(subset, plotFlag = False)



        # Plot
        if plotFlag:
            if eDim in subset.dims:
                # TODO: make this more robust. Currently assumes Euler dim is 1st dim.
                if subset[eDim].size > 1:
                    if pStyle is 'polar':
                        for EulerInd in range(0,subset[eDim].size):
#                                 print(f"*** Pol Geom: {subset[eDim][EulerInd].item()}")  # Need to pass all this through to plotters otherwise appears before plots!

                            # Set full title string here to pass to plotter.
                            # tString = f"{eDim}: {subset[eDim][EulerInd].item()}"
                            tString = f"Pol geom: {subset[eDim][EulerInd].item()}"

                            _ = sphSumPlotX(subset[EulerInd], pType = pType, backend = backend, facetDim = Etype, titleString = tString)

                    elif pStyle is 'grid':
                        print(f"Grid plot: {subset.attrs['jobLabel']}")

                        # Set data
                        subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)

                        if self.verbose['main']>2:
                            print(subset)

                        # If Phi is sliced, assumed 2D (Theta,E) plots
                        # NOTE - force real datatype here, otherwise get errors from complex residuals.
                        # TODO: why is this (after abs pType selection?) Should error check & clean up - this could cause issues sometimes...
                        if 'Phi' in selDims.keys():
                            # print(subset)
                            # subset.pipe(np.abs).plot(x='Theta', y=Etype, col=eDim, robust=True)
                            subset.plot(x='Theta', y=Etype, col=eDim, robust=True)
                        else:
                            # subset.pipe(np.abs).plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
                            subset.plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)

                        # plt.title(subset.attrs['jobLabel'])

            else:
                _ = sphSumPlotX(subset, pType = pType, backend = backend, facetDim = Etype)
