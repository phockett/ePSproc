"""
ePSproc plotting functions with Holoviews + Bokeh.

Aim: simple plotters for different datatypes, interactive in Jupyter Notebook + HTML formats.

15/07/20    Debugged, still pretty basic but running.
05/07/20    v1 in development.

See

 - https://epsproc.readthedocs.io/en/dev/tests/basicPlotting_dev_XC_030720.html
 - Plotting test notebooks (/tests/plottingDev) for more.
 - Dev code: http://localhost:8888/notebooks/github/ePSproc/epsproc/tests/plottingDev/basicPlotting_dev_280620.ipynb

Todo

- Plotting style mapping & options. Currently having HV issues here.
- Currently set only for XC datatypes from single dataSet, will want to enable stacking etc. here.
- Errorbar or spread plots, currently having issues getting these working for multidim data.

"""

import xarray as xr
# import holoviews as hv

# Optionals
# Additional plotters
# Seaborn for styles
try:
    import seaborn as sns
    snsFlag = True

except ImportError as e:
    if e.msg != "No module named 'seaborn'":
        raise
    print('* Seaborn not found, SNS styles not available. ')
    snsFlag = False

# Holoviews for plotting interface
try:
    import holoviews as hv
    from holoviews import opts
    hvFlag = True

except ImportError as e:
    if e.msg != "No module named 'holoviews'":
        raise
    print('* Holoviews not found, hvPlotters not available. ')
    hvFlag = False

# Set plotters & options.
def setPlotters(hvBackend = 'bokeh', width = 500):
    """
    Set some plot options - Seaborn style + HV defaults.

    May have some issues with scope here - TBC. Should just run on function import?
    UPDATE: now moved to module import.

    Parameters
    -----------
    hvBackend : str or list of strs, optional, default = 'bokeh'
        Backend(s) for holoviews to load.
        Can call bokeh, matplotlib or plotly

    width : int, optional, default = 500
        Setting for plot width, in pixels.

    """

    # Plotting libs
    # Optional - set seaborn for plot styling
    import seaborn as sns
    sns.set_context("paper")  # "paper", "talk", "poster", sets relative scale of elements
                            # https://seaborn.pydata.org/tutorial/aesthetics.html
    # sns.set(rc={'figure.figsize':(11.7,8.27)})  # Set figure size explicitly (inch)
                            # https://stackoverflow.com/questions/31594549/how-do-i-change-the-figure-size-for-a-seaborn-plot
                            # Wraps Matplotlib rcParams, https://matplotlib.org/tutorials/introductory/customizing.html
    sns.set(rc={'figure.dpi':(120)})

    from matplotlib import pyplot as plt  # For addtional plotting functionality
    # import bokeh

    # import holoviews as hv
    # from holoviews import opts

    # Set HV extension
    hv.extension(hvBackend)


    # Set global HV options
    # Setting frame_width here results in offset plots in layout - try setting later?
    # opts.defaults(opts.Curve(frame_width=500, tools=['hover'], show_grid=True, padding=0.01))
    opts.defaults(opts.Curve(width=width, tools=['hover'], show_grid=True, padding=0.01))

    # return hv.output.info()


# Convert "standard" XS dataarray to dataset format.
def hvdsConv(dataXS):
    """
    Basic conversion for XS data from Xarray to Holoviews.

    This will drop stacked Sym dims, and sum of Total to reduce - may not be appropriate in all cases?

    """

    ds = xr.Dataset({'sigma':dataXS.sel({'XC':'SIGMA'}).drop('XC'), 'beta':dataXS.sel({'XC':'BETA'}).drop('XC')})
    hv_ds = hv.Dataset(ds.unstack().sum('Total')) # OK - reduce Sym dims.

    return hv_ds, ds


# HV plotting routine for XS data
def XCplot(dataXS, lineDashList = {'L': 'dashed', 'M': 'solid', 'V': 'dashed'}):
    """
    Plot XC data using Holoviews.

    Currently optional stuff hard-coded here, will produce plots [sigma, beta] showing all data.
    Rather crude, needs some more style mapping.

    Parameters
    -----------
    dataXS : Xarray
        Xarray dataarray containing XC data in standard format.

    lineDashList : dict, optional, default = {'L': 'dashed', 'M': 'solid', 'V': 'dashed'}
        Set line types for calculation gauge.

    Returns
    --------
    layout : hv object


    Examples
    ---------

    >>> plotObj, _,_ = XCplot(dataXS[0])
    >>> plotObj

    """

    # Convert to HV dataset
    hv_ds, ds = hvdsConv(dataXS)

    # Set options
    # THIS NEEDS work!
#     lineDashList = {'L': 'dashed', 'M': 'solid', 'V': 'dashed'}
#     lineColorList = {'PU': 'blue', 'SU': 'red', 'All': 'green'}

    # This is working now... just need better cmapping
    dsPlotSet = hv.Layout()
    for vdim in ds.var():
        plotList = []
    #     for vdim in ds.var():
        for gauge in ds.Type:
    #         print(gauge.item())

            # Explicit looping here works for setting desired parameters independently
    #         dsPlotSet += hv_ds.select(Type=gauge.item()).to(hv.Curve, kdims=["Eke"], vdims=vdim, dynamic=False).opts(line_dash=lineList[gauge.item()]).overlay(['Cont'])
    #         plotList.append(hv_ds.select(Type=gauge.item()).to(hv.Curve, kdims=["Eke"], vdims=[vdim], dynamic=False).opts(line_dash=lineList[gauge.item()]).overlay(['Cont']))
            # Keep Type dim until *after* curve setting to allow for correct composition from list (otherwise will create Holomaps rather than curves)
                # With cmap also on type
        #         plotList.append(hv_ds.to(hv.Curve, kdims=["Eke"], vdims=[vdim, 'Type'], dynamic=False).select(Type=gauge.item()).opts(line_dash=lineDashList[gauge.item()], color=lineColorList[gauge.item()]).overlay(['Cont']))
                # Cmap on Cont
            plotList.append(hv_ds.to(hv.Curve, kdims=["Eke"], vdims=[vdim, 'Type'], dynamic=False).select(Type=gauge.item()).opts(line_dash=lineDashList[gauge.item()]).overlay(['Cont']))

        dsPlotSet += hv.Overlay(plotList)  #.groupby('Cont') #.collate()

#     (dsPlotSet + hv.Table(hv_ds)).cols(1)
    # return (dsPlotSet + hv.Table(hv_ds)).cols(1), hv_ds, ds
    # return (dsPlotSet).cols(1).opts(opts.Curve(frame_width=500)), hv.Table(hv_ds), hv_ds, ds
    return (dsPlotSet).cols(1), hv.Table(hv_ds), hv_ds, ds
