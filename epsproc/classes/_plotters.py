"""
ePSproc classes plot function wrappers

16/10/20    Methods for/from base class.

TODO: should be able to simplify with subselection logic removed (bit ad hoc currently due to organic growth of this!) and/or rewrite as decorators.


"""

from matplotlib import pyplot as plt  # For plot legends with Xarray plotter
import numpy as np # Needed only for np.nan at the moment
import xarray as xr

from epsproc import matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX, sphFromBLMPlot
from epsproc import lmPlot as lmPlotCore  # Hack rename here to prevent circular logic with local function - TODO: fix with core fn. reorg.
from epsproc.plot import hvPlotters
from epsproc.util.env import isnotebook
from epsproc.sphFuncs.sphConv import checkSphDims

# Import HV into local namespace if hvPlotters successful (can't access directly)
if hvPlotters.hvFlag:
    hv = hvPlotters.hv
    from epsproc.plot.util import showPlot

# ************** Plotters

def plotGetCro(self, pType = 'SIGMA', Erange = None, Etype = 'Eke', selDims = None, keys = None, backend = 'mpl'):
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

    selDims : str, optional, default = None
        Subselect dims, as a dictionary.
        E.g. to select a specific continuum symmetry set selDims = {'Cont':'Ag'}

    keys : list, optional, default = None
        Keys for datasets to plot.
        If None, all datasets will be plotted.

    backend : str, optional, default = 'mpl'
        Set plotter to use.

        - 'mpl' : Use Matplotlib/native Xarray plotter
        - 'hv' : use Holoviews via :py:func:`epsproc.plotters.hvPlotters.XCplot()`

    """
    # # Default to all datasets
    # if keys is None:
    #     keys = self.data.keys()
    keys = self._keysCheck(keys)

#         if self.jobs['jobStructure'] == 'subDirs':
    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
        # for m, item in enumerate(self.data[key]['XS']):

        # Set default to full range, same for all cases
        if Erange is None:
            Erange = [self.data[key]['XS'][Etype].min().data, self.data[key]['XS'][Etype].max().data]

        subset = self.data[key]['XS'].sel(selDims)  # SelDims here - uses dictionary as passed, OK for None too.
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
                    # plt.ylim([-1.5, 2.5])

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
def plotGetCroComp(self, pType='SIGMA', pGauge='L', pSym=('All','All'), Erange = None, Etype = 'Eke', Eshift = None, keys = None, backend = 'mpl', returnHandles = False):
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

    Eshift : int or float, optional, default = None
        Apply energy shift to results if set.

    keys : list, optional, default = None
        Keys for datasets to plot.
        If None, all datasets will be plotted.

    backend : str, optional, default = 'mpl'
        Set plotter to use.

        - 'mpl' : Use Matplotlib/native Xarray plotter
        - 'hv' : use Holoviews via :py:func:`epsproc.plotters.hvPlotters.XCplot()`

    returnHandles : bool, optional, default = False
        If true, return plot object and legend test list.

    NOTE: added backend options 27/10/20. CURRENTLY NOT WORKING for hv, due to data structure assumed in `hvPlotters.XCplot()`

    UPDATE 06/04/21: looking at this again, subselect on XC is one issue (can be fixed as per code plotGetCro()), but likely want to write this very differently for comparison case.

    """



    # Comparison plots over orbs

    # from matplotlib import pyplot as plt  # For legend
    lText = []

    # # Default to all datasets
    # if keys is None:
    #     keys = self.data.keys()
    keys = self._keysCheck(keys)

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

            lText.append(self.data[key]['jobNotes']['orbLabel'])

            # This is currently a bit dodgy, overwriting original data with in place update?
            # TODO: fix and test!
            if Eshift is not None:
                subset[Etype] += Eshift

            # Basic plot with Xarray
            # pltObj = subset.plot.line(x=Etype)
            # lText.append(self.data[key]['jobNotes']['orbLabel'])

            # Version with mpl or hv (crude - basically from XC plot above.)
            if backend == 'mpl':
                pltObj = subset.plot.line(x=Etype)
                # pltObj = subset.sel(XC=pType).plot.line(x=Etype)

                # if pType == 'SIGMA':
                #     plt.ylabel('XS/Mb')
                # else:
                #     plt.ylabel(r"$\beta_{LM}$")

            # NOTE: added backend options 27/10/20. CURRENTLY NOT WORKING for hv, due to data structure assumed in `hvPlotters.XCplot()`
            # TODO: change formatting to use XCplot and/or set up new subfunction.
            if backend == 'hv':
                pltObj, *_ = hvPlotters.XCplot(subset, kdims = Etype)

                print(f"\n*** {jobLabel}")
                # Check if notebook and output
                if isnotebook():
                    display(pltObj) # Use IPython display(), works nicely for notebook output
                    # Q: is display() always loaded automatically? I think so.
                else:
                    pltObj  # Not yet tested - not sure what the options are here for hv/Bokeh.
                            # May just want to return this object (hv layout)?


        # Label with orb_sym
#                 lText.append(self.dataSets[key]['XS'][m].attrs['fileBase'].rsplit('/',maxsplit=1)[0])

        # lText.append(f"Orb {self.dataSets[key]['XS'][m].attrs['orbInfo']['orbN']} ({self.dataSets[key]['XS'][m].attrs['orbInfo']['orbSym'][0]})")

        # lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])


    # Update legend etc.
    if backend == 'mpl':
        plt.legend(lText)

        if pType == 'SIGMA':
            plt.ylabel('XS/Mb')
        else:
            plt.ylabel(r"$\beta_{LM}$")
            # plt.ylim([-1.5, 2.5])

    if returnHandles:
        return pltObj, lText




def lmPlot(self, Erange = None, Etype = 'Eke', dataType = 'matE', xDim = None, keys = None, refDataKey = None, reindexTol = 0.5, reindexFill = np.nan, setPD = True, **kwargs):
    """
    Wrapper for :py:func:`epsproc.lmPlot` for multijob class. Runs lmPlot() for each dataset.

    Parameters
    ----------
    Erange : list of int or float, optional, default = None
        Set plot range [Emin, Emax]. Defaults to full data range if not set.

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

    # # Default to all datasets
    # if keys is None:
    #     keys = list(self.data.keys())
    keys = self._keysCheck(keys)

    # Set lmPlotOpts
    # Check passed args vs. self.lmPlotOpts and overwrite
    if kwargs:
        for key, value in kwargs.items():
            self.lmPlotOpts[key] = value

    # Check Etype exists, it may not for some data types.
    # if not Etype in self.data[keys[0]].keys():
    if not Etype in list(self.data[keys[0]][dataType].coords):
        Etype = None

    # Set default to full range of 1st dataset, keep same for all cases
    # TODO: set per dataset?
    if (Erange is None) and Etype:
        Erange = [self.data[keys[0]][dataType][Etype].min().data, self.data[keys[0]][dataType][Etype].max().data]

    # Set ref dataset if required
    if refDataKey is not None:
        refData = self.data[refDataKey][dataType]

        if Etype == 'Ehv':
            refData = refData.swap_dims({'Eke':'Ehv'})

        if Etype:
            refData = refData.sel(**{Etype:slice(Erange[0], Erange[1])})  # Case for slicing on Etype

        refData.attrs = self.data[refDataKey][dataType].attrs # Propagate atrrs.

    else:
        refData = None


    # Loop over datasets
    for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works

        # Init empty list for daPlotpd data
        # 21/10/20 - should only be single item per key in current data schema, so just set directly below
        # if setPD:
        #     self.data[key]['daPlotpd'] = []

        # for m, item in enumerate(self.data[key]['matE']):

        # More elegant way to swap on dims?
        if Etype == 'Ehv':
        # Subset before plot to avoid errors on empty array selection!
            subset = self.data[key][dataType].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        elif Etype == 'Eke':
            subset = self.data[key][dataType].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

        else:
            subset = self.data[key][dataType]  # Case for no Etype dim.

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
        # 21/10/20 - should only be single item per key in current data schema, so just set directly
        if setPD:
            # self.data[key]['daPlotpd'].append(daPlotpd)  # Set to include None cases to keep indexing. Should set as dict instead?
            self.data[key]['daPlotpd'] = daPlotpd

def ADMplot(self, dataType = 'ADM', xDim = 't', Etype='t', col = None, **kwargs):
    """
    Wrap BLMplot() for ADMs. Thin wrapper with some ADM-specific defaults.

    TODO: make this good.
    """
    # Pass using locals()
    # self.BLMplot(**locals())  # Almost neat, but needs logic to remove self and unpack kwargs

    # Pass explicitly
    self.BLMplot(dataType = dataType, xDim = xDim, Etype = Etype, col = col, **kwargs)


def BLMplot(self, Erange = None, Etype = 'Eke', dataType = 'AFBLM',
            xDim = None, selDims = None, col = 'Labels', row = None,
            thres = None, keys = None, verbose = None,
            backend = 'xr', overlay = None,
            **kwargs):
    """
    Basic BLM line plots using Xarray plotter.

    See https://epsproc.readthedocs.io/en/latest/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html

    Similar to :py:func:`epsproc.BLMplot`, may change to simple wrapper, but some differenences in terms of dim handling here.

    For more flexibility, use `self.lmPlot`.

    TODO: update BLMplot to support more datatypes, and implement here instead.

    TODO: fix dim handling and subselection, see old plotting code.

    24/11/21: quick additions, override printing with "verbose", and added backend option for XR or Holoviews plotters.
                Note this currently uses hvplot functionality, see https://hvplot.holoviz.org/user_guide/Gridded_Data.html.
                UPDATE: currently not working due to unhandled dims at Holomap stack - see tmo-dev for method.
    05/06/21: added **kwargs pass to Xarray line plot
    03/02/21: added col, row arguments for flexibility on calling. Still needs automated dim handling.

    """

    if backend == 'hv':
        # self._hvBLMplot(**locals())
        self._hvBLMplot(**{k:v for k,v in locals().items() if k not in ['self','kwargs']}, **kwargs)

    # Set xDim if not passed.
    if xDim is None:
        xDim = Etype

    if verbose is None:
        verbose = self.verbose

    # # Default to all datasets
    # if keys is None:
    #     keys = list(self.data.keys())
    keys = self._keysCheck(keys)

    # Set default to full range of 1st dataset, keep same for all cases
    # TODO: set per dataset? Now in subset fn.
    # if Erange is None:
    #     Erange = [self.data[keys[0]][dataType][Etype].min().data, self.data[keys[0]][dataType][Etype].max().data]

    # Loop over datasets
    plotList = {}  # For HV plot stacking
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

        # Threshold results. Set to check along Eke dim, but may want to pass this as an option.
        # Set squeeze to True also
        subset = matEleSelector(subset, thres=thres, inds = selDims, dims = Etype, sq = True)

        if subset.any():
            # THIS IS SHIT
            # if hasattr(subset, 'XSrescaled'):
            # #     print(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, XS")
            #     # list(set(subset.dims) - {row} - {col})
            #     rowXS = row
            #     colXS = col
            #     if 'BLM' in row:
            #         rowXS = None
            #     if 'BLM' in col:
            #         colXS = None
            #     subset.XSrescaled.real.plot(x=Etype, col=colXS, row=rowXS)   # UGH THESE DIMENSION ARE NOT THE SAME OF COURSE SO WILL POTENTIALLY BREAK. SHITTY CODE AGAIN.
            #     # plt.title(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, XS")

            if verbose:
                try:
                    print(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, {dataType}")
                except KeyError:
                    print(f"Dataset: {key}, {dataType}")

            if backend == 'xr':
                # TODO: add some logic here, sort or switch on flag or number of dims?
                # subset.real.plot(x=Etype, col='Labels', row='BLM')  # Nice... should give line plots or surfaces depending on dims
                # subset.where(subset.l>0).real.plot.line(x=Etype, col=col, row=row)  # Set to line plot here to stack BLMs, also force l>0 - THIS SUCKS, get an empty B00 panel.
                subset.real.plot.line(x=Etype, col=col, row=row, **kwargs)
                # plt.title(f"Dataset: {key}, {self.data[key]['jobNotes']['orbLabel']}, {dataType}")

            # VERY ROUGH, based on not great tmo-dev code.
            # Currently no explicit dim handling, and auto-rename to ensure working Holomap output.
            # TODO: update/replace with routines in basicPlotters.BLMplot(), see also hvPlotters.curvePlot()
            elif backend == 'hv':
                try:
                    # plotList[key] = (subset.real.hvplot.line(x=Etype, col=col, row=row, **kwargs))
                    hvDS = hv.Dataset(subset.real.rename(str(key)))  # Convert to hv.Dataset, may need to rename too.
                    # hvDS = hv.Dataset(subset.real.rename('AFBLM'))
                    plotList[key] = hvDS.to(hv.Curve, kdims=Etype)   # TODO: add dim handling, will just autostack currently.

                    # plotList[key] = hvPlotters.curvePlot(dataXR, kdims = Etype, returnPlot = True, renderPlot = False, **kwargs) # Curve plot version - needs better dim handling too!

                except NotImplementedError as e:
                    # if e.msg == "isna is not defined for MultiIndex":
                    if e.args[0] == "isna is not defined for MultiIndex":   # Message as first arg - may be version specific?
                        print("Unstacking data for plotting")
                        # plotList[key] = (subset.real.unstack().hvplot.line(x=Etype, col=col, row=row, **kwargs))   # Should be OK for individual plots, but can't pass to holomap?
                        hvDS = hv.Dataset(subset.real.unstack().rename(str(key)))  # Convert to hv.Dataset, may need to rename too.
                        # hvDS = hv.Dataset(subset.real.unstack().rename('AFBLM'))  # Convert to hv.Dataset, may need to rename too.
                        plotList[key] = hvDS.to(hv.Curve, kdims=Etype)   # TODO: add dim handling, will just autostack currently.
                    else:
                        print(f"Caught unhandelable exception {e}")


    if backend == 'hv':
        # 31/10/22 - debugging issues with HoloMap - see ep.basicPlotters.BLMplot(), which just uses overlays?
        try:
            # Code from showPlot()
            hvMap = hv.HoloMap(plotList)  # May need to handle dims here?
            if self.__notebook__ and (backend == 'hv'):
                if overlay is None:
                    display(hvMap)  # If notebook, use display to push plot.
                else:
                    display(hvMap.overlay(overlay))
        except:
            print("Failed to compose HoloMap, returning plot dictionary.")
            return plotList


# New BLMplot for hv backend only, but may end up more general
# def _hvBLMplot(self, **kwargs):  # This works, but all args nested in kwargs dict
def _hvBLMplot(self, Erange = None, Etype = 'Eke', dataType = 'AFBLM',
            xDim = None, selDims = None, col = 'Labels', row = None,
            thres = None, keys = None, verbose = None,
            backend = 'xr', overlay = None,
            **kwargs):
    """
    Plot BLM data with Holoviews. Subfunction for self.BLMplot, but also aim to implement better and general dim handling.

    Method follows dev code at https://phockett.github.io/ePSdata/OCS-preliminary/OCS_orbs8-11_AFBLMs_VM-ADMs_140122-JAKE_tidy.html#AFBLMs-for-aligned-case
    Subset data to Xarray Dataset then plot, rather than compose multiple plots and stack as per old method - should be simpler and more flexible & robust.

    TODO: clarify arg passing, currently just dump locals to kwargs from main self.BLMplot call.

    """

    print(locals())
    #*** Var checks - currently as per self.BLMplot
    # Set xDim if not passed.
    if xDim is None:
        xDim = Etype

    if verbose is None:
        verbose = self.verbose

    if not ('pType' in locals()) or (pType is None):
        pType = 'r'

    if not ('renderPlot' in locals()):
        renderPlot = True

    if not ('returnPlot' in locals()):
        returnPlot = True

    # # Default to all datasets
    # if keys is None:
    #     keys = list(self.data.keys())
    keys = self._keysCheck(keys)


    # UPDATED ROUTINE FROM https://phockett.github.io/ePSdata/OCS-preliminary/OCS_orbs8-11_AFBLMs_VM-ADMs_140122-JAKE_tidy.html#AFBLMs-for-aligned-case
    # STACK TO XR FOR ALL DIM HANDLING, then plot!!!
    # import xarray as xr
    pDict = []
    for key in keys:
    #     subset = ep.matEleSelector(data.data[key][dataType], thres=thres, inds = selDims, dims = [Etype, 't'], sq = True)
        subset = matEleSelector(self.data[key][dataType], thres=thres, inds = selDims, dims = Etype, sq = True)
        subset.name = 'BLM'

        subset = plotTypeSelector(subset, pType = pType)

    #     pDict.append(subset.expand_dims({'Orb':[key]}))
    #     pDict.append(subset.expand_dims({'Orb':[data.data[key][dataType].attrs['jobLabel']]}))
        pDict.append(subset.expand_dims({'Orb':[self.data[key]['jobNotes']['orbLabel']]}))


    xrDS = xr.concat(pDict, 'Orb')

    #*** Also check BLM dims for stacking?
    # Note this sets xrDS.attrs.harmonics, and these must all be the same.
    checkSphDims(xrDS)

    #*** Holoviews init
    hvDS = hvPlotters.hv.Dataset(xrDS.unstack(xrDS.attrs['harmonics']['stackDim']))  # Need to unstack, otherwise sometimes get empty plots?
                                                                                      # Seems to be issue with tuple/MultiIndex case.

    hvObj = hvDS.to(hvPlotters.hv.Curve, kdims=Etype)  # .overlay(xrDS.attrs['harmonics']['dimList'])   # OK

    if renderPlot:
        return showPlot(hvObj.overlay(xrDS.attrs['harmonics']['dimList']), returnPlot = returnPlot, __notebook__ = isnotebook())  # Currently need to pass __notebook__?
    else:
        return hvObj


# # Plot PADs from mat elements (MF) or BLMs
# def padPlot(self, selDims = {}, sumDims = {}, Erange = None, Etype = 'Eke',
#                 keys = None, dataType = 'TX',
#                 pType = 'a', pStyle = 'polar', plotFlagGlobal = True,
#                 backend = 'mpl'):
#
#     """
#     Plot PADs from mat elements (MF) or BLMs.
#     """
#
#     # # Default to all datasets
#     # if keys is None:
#     #     keys = self.data.keys()
#     # else:
#     #     if not isinstance(keys, list):   # Force list if single item passed
#     #         keys = [keys]
#     keys = self._keysCheck(keys)
#
#     # Loop over datasets
#     for key in keys:
#         # for m, item in enumerate(self.dataSets[key]['TX']):
#         plotFlag = True
#
#         # # More elegant way to swap on dims?
#         # subset = self.data[key][dataType]
#         #
#         # if Etype == 'Ehv':
#         # # Subset before plot to avoid errors on empty array selection!
#         #     subset = subset.swap_dims({'Eke':'Ehv'})
#         #
#         # # Slice on Erange with/without step size - bit ugly.
#         # if Erange is not None:
#         #     if len(Erange) == 3:
#         #         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword
#         #     else:
#         #         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})
#         subset = self.Esubset(key = key, dataType = dataType,  Erange = Erange, Etype = Etype)
#
#         # Handling for Euler or Labels dim
#         eDim = 'Euler'
#         if hasattr(subset, 'Labels'):
# #                 if 'Labels' in selDims:
#             if eDim in subset.dims:
#                 subset = subset.swap_dims({'Euler':'Labels'})
#
#             eDim = 'Labels'
#
#         # if selDims is not None:  # Now set as empty dict to avoid issues later.
#         # TODO: smarter dim handling here - selDims may not be consistent over all dataSets!
#         if selDims:
#             subset = subset.sel(selDims)
#
#
#         # Compute PADs if not already present in dataset
#         # TODO: more careful type checking here, should have general dist case?
#         # 16/01/21 moved up to allow for dim checking for grid case - BUT THIS IS SHIT CODE, LOTs of this already in sph plot fns.
#         if dataType not in ["TX", "wigner"]:
#             subset, _ = sphFromBLMPlot(subset, plotFlag = False)
#             subset = subset.sum('LM')
#
#
#         if self.verbose['main']>2:
#             print(subset.dims)
#
#         # Check dimensionality - sum over Sym & it if necessary
#         if (subset.ndim > 3):  # Need some additional dim logic here!
#             print(f"Found dims {subset.dims}, summing to reduce for plot. Pass selDims to avoid.")
#             if 'Sym' in subset.dims:
#                 subset = subset.sum('Sym').squeeze()
#             if 'it' in subset.dims:
#                 subset = subset.sum('it').squeeze()
#
#             # Check result is OK
#             # TODO: check for specific dims here
#             if (subset.ndim > 3): # and not (eDim in subset.dims):
#                 # Try squeezing in case of singleton dims
#                 subset = subset.squeeze()
#
#                 # Define non-parseable case - should make this more general.
#                 # This also isn't correct for handling BLM vs. gridded data types.
#                 testDims = list(set(subset.dims) - set(eDim))
#                 if self.verbose['main']>2:
#                     print(testDims)
#
#                 if (len(testDims) > 3):
#                     print(f"*** ERROR: dataset {self.data[key][dataType].jobLabel}  ndims = {subset.ndim}, dims = {subset.dims}. Skipping MFPAD plotting.")
#                     plotFlag = False
#
#             # else:
#             #     pass
#
#         # Propagate attrs
#         subset.attrs = self.data[key][dataType].attrs
#         #
#         # # Compute PADs if not already present in dataset
#         # # TODO: more careful type checking here, should have general dist case?
#         # if dataType not in ["TX", "wigner"]:
#         #     subset, _ = sphFromBLMPlot(subset, plotFlag = False)
#
#
#         # Plot
#         if plotFlag and plotFlagGlobal:
#             if eDim in subset.dims:
#                 # TODO: make this more robust. Currently assumes Euler dim is 1st dim.
#                 if subset[eDim].size > 1:
#                     if pStyle is 'polar':
#                         for EulerInd in range(0,subset[eDim].size):
# #                                 print(f"*** Pol Geom: {subset[eDim][EulerInd].item()}")  # Need to pass all this through to plotters otherwise appears before plots!
#
#                             # Set full title string here to pass to plotter.
#                             # tString = f"{eDim}: {subset[eDim][EulerInd].item()}"
#                             tString = f"Pol geom: {subset[eDim][EulerInd].item()}, ploType: {pType}"
#
#                             _ = sphSumPlotX(subset[EulerInd], pType = pType, backend = backend, facetDim = Etype, titleString = tString)
#
#                     elif pStyle is 'grid':
#                         print(f"Grid plot: {subset.attrs['jobLabel']}, dataType: {dataType}, plotType: {pType}")
#
#                         # Set data
#                         subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)
#
#                         if self.verbose['main']>2:
#                             print(subset)
#
#                         # If Phi is sliced, assumed 2D (Theta,E) plots
#                         # NOTE - force real datatype here, otherwise get errors from complex residuals.
#                         # TODO: why is this (after abs pType selection?) Should error check & clean up - this could cause issues sometimes...
#                         if 'Phi' in selDims.keys():
#                             # print(subset)
#                             # subset.pipe(np.abs).plot(x='Theta', y=Etype, col=eDim, robust=True)
#                             subset.plot(x='Theta', y=Etype, col=eDim, robust=True)
#                         else:
#                             # subset.pipe(np.abs).plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
#                             subset.plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
#
#                         # plt.title(subset.attrs['jobLabel'])
#
#             else:
#                 # 15/01/21: developing/debugging here - this case seems to be called even with data WITH EULERs?  Why?
#                 # tString = f"Pol geom: {subset[eDim][EulerInd].item()}, ploType: {pType}"
#                 tString = 'No Euler test'
#
#                 if pStyle is 'polar':
#                     _ = sphSumPlotX(subset, pType = pType, backend = backend, facetDim = Etype, titleString = tString)
#
#                 elif pStyle is 'grid':
#                     # Autoset dims from non (Theta,Phi) Dims
#                     rowDim = list(set(subset.dims) - set(('Theta','Phi',Etype)))
#
#                     # Set data
#                     subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)
#
#                     if len(rowDim) > 1:
#                         print(f"*** ERROR: gridplot skipped for  {self.data[key][dataType].jobLabel}, rowDims = {rowDim} too many dims.")
#
#                     elif len(rowDim) == 1:
#                         print(f"Grid plot: {self.data[key][dataType].attrs['jobLabel']}, dataType: {dataType}, plotType: {pType}, dims ({rowDim}, {Etype})")
#                         subset.plot(y='Theta',x='Phi', row=rowDim, col=Etype, robust=True)
#
#                     else:
#                         print(f"Grid plot: {self.data[key][dataType].attrs['jobLabel']}, dataType: {dataType}, plotType: {pType}, dims ({Etype})")
#                         subset.plot(y='Theta',x='Phi', col=Etype, robust=True)
#
#     # Return data for further use
#     if not plotFlagGlobal:
#         return subset

# padPlot v2, rewritten 19/01/21
# This has better dim handling, and simpler logic, but still needs some work.
def padPlot(self, selDims = {}, sumDims = {'Sym','it'}, Erange = None, Etype = 'Eke',
                keys = None, dataType = 'TX', facetDims = None, squeeze = False, reducePhi = None,
                pType = None, pStyle = 'polar', returnFlag = False, plotDict = 'plots',
                backend = 'mpl', plotFlag = True):

    """
    Plot I(theta,phi) data from BLMs or gridded datasets.

    reducePhi : optional, default = None
    This allow phi selection or summation for parameters which required Ylm expansion before plotting.
        Pass 'sum' to sum over phi before plotting.
        Pass a value to select.


    plotFlag : bool, optional, default = True
        Set plotFlag=False bypass for plotter object return only.

    pType : str, optional, default = None
        Plot type, defaults to 'a' (abs) except for dataType = 'TX', 'a2'.
        See :py:func:`epsproc.sphPlot.plotTypeSelector` for options.

    TODO: fix dim handling for pl case, need to pass facetDim != None.
    TODO: return plot objects. Probably to self.data[key][pStyle], or as dictionary of plots per run with data? (Cf. PEMtk plotters.)
        23/04/22: added plot data and object returns to self.data[key][plotDict][pStyle]

    """


    # Default to all datasets
    keys = self._keysCheck(keys)

    if pType is None:
        pType = 'a'

        if dataType == 'TX':
            pType = 'a2'

    # Default facetDims, (Eulers, Eke)
    # Should test these? YES - it's done later, per key.
    if facetDims is None:
        facetDims = ['Labels', Etype]

    if len(facetDims) == 1:
        facetDims.append(None)
    if len(facetDims) > 2:
        print(f"Too many facetDims (max=2), using only {facetDims[0:2]}.")


    # Loop over datasets
    for key in keys:
        # for m, item in enumerate(self.dataSets[key]['TX']):
        # plotFlag = True

        # Set data & slice on E
        subset = self.Esubset(key = key, dataType = dataType,  Erange = Erange, Etype = Etype)

        # Handling for Euler or Labels dim - UGLY
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

        # If theta and/or phi present, assume angular gridded data
        # If not present, assumed LM parameters, and expand on grid
        if not (('Theta' in subset.dims) or ('Phi' in subset.dims)):
            subset, _ = sphFromBLMPlot(subset, plotFlag = False)
            subset = subset.sum('LM')

        if sumDims:
            sumDimsCheck = set(subset.dims)&{*sumDims}  # This checks sumDims are present, otherwise will throw an error.
            print(f"Summing over dims: {sumDimsCheck}")
            subset = subset.sum(sumDimsCheck)

        if reducePhi:
            if reducePhi == 'sum':
                subset = subset.sum('Phi')
            else:
                subset = subset.sel({'Phi':reducePhi})

        if squeeze:
            subset = subset.squeeze()


        #***** Check dims
#         # Handling for Euler or Labels dim - UGLY... now moved above to allow for selDims in Xarray > 0.15 (treats non-dimensional coords differently)
#         eDim = 'Euler'
#         if hasattr(subset, 'Labels'):
# #                 if 'Labels' in selDims:
#             if eDim in subset.dims:
#                 subset = subset.swap_dims({'Euler':'Labels'})
#
#             eDim = 'Labels'

        #
#         plotDims = set(subset.dims) - set(facetDims)

        # Check facetDims exist, otherwise may get groupby errors (may be cleaner to use try/except here?)
        # NOTE this will still produce errors in some cases (0 dims matched)
        facetDimsCheck = list(set(subset.dims)&{*facetDims})
        groupFlag = True
        if len(facetDimsCheck) == 0:
            print(f'***Error: missing dims {facetDims}.')
        elif len(facetDimsCheck) == 1:
            groupFlag = False
            facetDimsCheck.append(None)

        extraDims = set(subset.dims) - {*facetDimsCheck,*sumDims} - {'Theta','Phi'}  # Check for outstanding dims, this will return an empty set if all dims accounted for here

        if extraDims:
            print(f"Found additional dims {extraDims}, summing to reduce for plot. Pass selDims to avoid.")

            for dim in extraDims:
                subset = subset.sum(dim)  #.squeeze()

        # Ensure attrs propagated
        subset.attrs = self.data[key][dataType].attrs

        print(f"Plotting from self.data[{key}][{dataType}], facetDims={facetDimsCheck}, pType={pType} with backend={backend}.")

        # GROUPBY AND PLOT?  NOT SURE HOW BEST TO HANDLE MULTIPLE DIMS HERE? Can pass 1 facet dim to SPH plotter already.
        # TODO: decide MAX 4D. Then reduce > groupby > facet to existing plotter.
        # TODO: general way to handle more dims?
        #
        # 12/08/22: added group by-pass for single facet dim case.
        #           Seems to be OK for Plotly case, but NOT WORKING CORRECTLY FOR MPL CASE HOWEVER, still get single fig per item?

        if pStyle == 'polar':
            if groupFlag:
                for groupLabel, item in subset.groupby(facetDimsCheck[0]):
    #                 tString = f"Pol geom: {item.item()}, plotType: {pType}"
                    tString = f"{facetDimsCheck[0]}: {groupLabel}, plotType: {pType}"
                    pltObj = sphSumPlotX(item, pType = pType, backend = backend, facetDim = facetDimsCheck[1], titleString = tString, plotFlag=plotFlag, verbose = self.verbose['sub'])
                    # _ = sphSumPlotX(item, pType = pType, backend = backend, facetDim = facetDimsCheck[1], titleString = tString)
            else:
                # For single dim case just pass to existing plotter for subplots,
                tString = f"plotType: {pType}"
                pltObj = sphSumPlotX(subset, pType = pType, backend = backend, facetDim = facetDimsCheck[0], titleString = tString, plotFlag=plotFlag, verbose = self.verbose['sub'])


        elif pStyle == 'grid':
            if 'jobLabel' in subset.attrs.keys():
                label = subset.attrs['jobLabel']
            else:
                label = (key,dataType)

            print(f"Grid plot: {label}, dataType: {dataType}, plotType: {pType}")

            # Set data
            subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)

            try:
                if reducePhi:
                    # subset.plot(x='Theta', y=Etype, col=eDim, robust=True)
                    # subset.plot(x='Theta', y=facetDimsCheck[0], col=facetDimsCheck[1], robust=True)  # This might fail for only one facetDim
                                                                                                        # This did work initially, but not later - dim ordering always goes to alphabetical with set selection above
                    pltObj = subset.plot(x='Theta', y=Etype, col=list({*facetDimsCheck}-{Etype})[0], robust=True)  # Force E dim to y
                else:
                    pltObj = subset.plot(y='Theta',x='Phi', row=facetDimsCheck[1], col=facetDimsCheck[0], robust=True)  # This might fail for only one facetDim. DIM ORDERING NOT PRESERVED

            except ValueError as e:
                print(f"*** Error {e} for grid plotter with facetDims={facetDimsCheck}: try passing selDims, sumDims or facetDims manually.")  # Gives "ValueError: IndexVariable objects must be 1-dimensional" for singleton facet dim case (gets squeezed out here?)


        # Return data?
        if returnFlag:
            if not plotDict in self.data[key].keys():
                self.data[key][plotDict] = {}

            self.data[key][plotDict][dataType] = {'pData':subset,
                                                   pStyle:pltObj}

            if self.verbose:
                print(f"Set plot to self.data['{key}']['{plotDict}']['{dataType}']['{pStyle}']")

    # # Return data? Note this is only set for final key value at the moment.
    # if returnFlag:
    #     return subset
