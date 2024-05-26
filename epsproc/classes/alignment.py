"""
Basic class for handling alignment (ADM) data

23/03/24   v1, basic implementation using ePSmultiJob as base class.
               Implements file scanner, basic handling and normalisation, basic plotting (hv) and basic subset routine.
               Code based on previous manual routines, e.g. OCS dev work, https://phockett.github.io/ePSdata/OCS-preliminary/OCS_orbs8-11_AFBLMs_VM-ADMs_140122-JAKE_tidy-replot-200722_v5.html

For initial testing, see notebook ePS/N2O/epsman2024/proc/N2O_ADMs_test_200324.ipynb


"""


from epsproc.classes.base import ePSbase
from epsproc.classes.multiJob import ePSmultiJob
from epsproc.util.conversion import datasetStack
from epsproc.IO import getFiles
from epsproc.sphCalc import setADMs

import pandas as pd
import numpy as np
from scipy.io import loadmat
from pathlib import Path

# Import HV into local namespace if hvPlotters successful (can't access directly)
from epsproc.plot import hvPlotters
if hvPlotters.hvFlag:
    hv = hvPlotters.hv
    from epsproc.plot.util import showPlot

class ADM(ePSmultiJob):
    """
    Testing class-based ADM IO and methods.

    v1 20/03/24
    """

    pass

    # Basic dir + subdir scanner
    # NOTE ep.getFiles has no subdir support!

    def scanDirs(self):

        # Quick subdir scan
        super().scanDirs()

        # Now just run getFiles per subdir...
        fileList = []
        fileDict = {}
        for item in self.jobs['jobDirs']:
            # Basic list
            subFiles = getFiles(fileBase =  item, fType=self.jobs['ext'])
            fileList.append(subFiles)
            # Set dict by item names...
            fileDict[item.name] = subFiles


        self.fileList = fileList
        self.fileDict = fileDict


    def loadData(self, keyType = 'fileName',
                 fReader=None, fType = 'text',
                 normType = None, renorm = False,
                 addPop = True, addS = True):
        """
        Load ADM data from file(s).

        NOTE: For `Tunstack` case currently assumes dir and file structure as:
            root dir - AKQS dir, e.g. 'A200' - files per T, e.g. 'A200_10K.txt'
            Variables are parsed according to this scheme (for 2024 N2O data).

        For older datasets (single dir/alternative naming schema) may need converters or regex methods.
            Using `keyType='fileName'` should generally read data, but may not pull labels correctly currently, and may fail to convert to Xarray.

        Time points are searched for in PARENT dir.


        Parameters
        ----------

        keyType : str, optional, default = 'fileName'
            Set data type and load style (needs work).
            Currently:
            - 'fileName', load and set to dictionary with filename keys.
                (Note this will set self.data[dir] and stack ADMs by file.)
            - 'T', load and set to dictionary with T (temperature) keys.
            - 'Tunstack', load and set to dictionary with T (temperature) keys, alternative formatting.
            - 'setADMsT', load and set to XR dataArray with T (temperature) keys.

        fReader : function, optional, default = None
            Pass file read function.
            If None, defaults will be set by file ext.
            - For .txt, .csv, use pd.read_csv
            - For .mat, .dat, use scipy.io.loadmat

        fType : str, optional, default = 'text'
            Used to explicitly set file type if fReader set.
            Otherwise set automatically by file extension.
            Note this only controls some extra data formatting options.

        normType : str, optional, default = None
            Set normType (via self.normDict).
            If None, skip this, and also ignore addPop.

        renorm : bool, optional, default = True
            Apply additional renormalisation by (2*K+1)/8*pi^2

        addPop : bool, optional, default = True
            Add K=0 population term?
            Only applied if normType is set.

        addS : bool, optional, default = True
            Add S=0 terms?
            Default=False.

        TODO:
        - More general searching and file handling.
            - Not all methods use norm or renorm.
            - 15/05/24 added general handling for multiple files per dir.
                Should test further, see other libs for alternative data structures...?
                Currently stack per dir, files labelled as dict items or in xr.Dataset.
        - Automate addS, currently setADMs throws errors if not correct dims set.

        """

        # Set file reader by type
        if fReader is None:
            if self.jobs['ext'] in ['.txt','.csv']:
                fReader = pd.read_csv
                fType = 'text'
            elif self.jobs['ext'] in ['.mat','.dat']:
                fReader = loadmat
                fType = 'matlab'
            else:
                print(f"*** File type {self.jobs['ext']} not recognised.")
                print("Skipping file reading. To force, pass fReader=<file read function>. E.g. `fReader=pd.read_csv` for default case.")
                return None

        # Set renorm/population term
        # If None this will NOT be set
        self.normDict(normType = normType)

        # Try loading t-data from root dir
        # Assume only one file?
        # NOTE - currently only used for setADMs style IO method.
        tFiles = list(self.job['fileBase'].parent.glob('t*'))
        if tFiles:
            print(f"Found t-index {tFiles[0]}")
            tIndex = pd.read_csv(tFiles[0], header=None).to_numpy().squeeze()
#             print(tIndex)
        else:
            tIndex = None

        Tcols = []  # Set this for output, not used in all cases.

        # Scan data files with fReader - loop over subdirs and files
        # TODO: add file name as XR attrs.
        if keyType == 'fileName':
            dataDict = {k:{Path(f).name:fReader(f) for f in item} for k,item in self.fileDict.items()}

            # Try assigning ADMs per dir...
            for k,item in dataDict.items():
                xrData = {}
                if fType == 'matlab':
                    # For Matlab case, have items for things in single-file case
                    # xrData = {k2:setADMs(ADMs = item2['ADM'], t=item2['time'].squeeze(),
                    #             KQSLabels = item2['ADMlist'], addS = addS) for k2,item2 in item.items() if 'time' in item2.keys()}

                    for k2,item2 in item.items():
                        if 'time' in item2.keys():
                            xrData[k2] = {'ADM':setADMs(ADMs = item2['ADM'], t=item2['time'].squeeze(),
                                                KQSLabels = item2['ADMlist'], addS = addS)  }

                    # Set also to standard self.data[key][dataType] style
                    # This is used by existing core functionality
                    # TODO: fix to use something like xrDA, xrDS, dataDict = datasetStack(ADMin.data['alignment']['ADM'], dataType='ADM', stackDim = 'file', keys = ADMin.data['alignment']['ADM'].keys())
                    #       Currently this fails, since the ordering is incorrect - should reorder to data[dir][fileName]['ADM'] for this.
                    # self.data[k] = xrData

                    # Stack to single DA/DS per dir...
                    # Use existing function, but note need to set coords curently
                    fKeys = list(xrData.keys())
                    fNames = [item.rstrip(self.jobs['ext']) for item in fKeys]
                    xrDA, xrDS, dataDict = datasetStack(xrData, dataType='ADM', stackDim = 'file', keys = xrData.keys())

                    # TODO - store XR per dir too...?
                    self.xr = xrDA.assign_coords({"file":fNames})
                    self.xrDS = xrDS
                    # self.xr.Temp.attrs['units'] = 'K'
                    # self.xr = self.xr.sortby('Temp')

                    # Set data with dir as key
                    # self.data[k] = {'ADM':self.xr}
                    self.data[k] = {'ADM':xrDA.assign_coords({"file":fNames})}  # Set again here as above will be overriden!

                else:
                    print(f"No data formatting implemented for fType '{fType}' from dir '{k}'. Please run ep.setADMs manually to reformat")

        if keyType == 'T':
            dataDict = {k:{f.split('_')[-1].strip(self.jobs['ext']):pd.read_csv(f) for f in item} for k,item in self.fileDict.items()}

        # Case for dirs per Temp.
        # CURRENTLY ONLY SUPPORTS pd.read_csv.
        if keyType == 'Tunstack':
            # Try better formatting for T data...
            # See also https://stackoverflow.com/a/21232849 for ideas/cleaner methods
            # Xr directly?
            # See https://docs.xarray.dev/en/stable/generated/xarray.Dataset.from_dict.html

            dataDict = {}

            for k,item in self.fileDict.items():
#                 pd.read_csv
#                 Temp = f.split('_')[-1].strip(self.jobs['ext']
#                 Tarray = np.array([pd.read_csv(f).to_numpy() for f in item])  # Numpy array

                Tframe = pd.concat([pd.read_csv(f, names=[f.split('_')[-1].strip(self.jobs['ext'])], header=None) for f in item], axis=1) #, ignore_index=True)
#                 Tframe = pd.concat([pd.read_csv(f) for f in item], axis=1, ignore_index=True)
                     # Adding .rename(columns=[f.split('_')[-1].strip(self.jobs['ext'])]) fails.

                Tcols = [f.split('_')[-1].strip(self.jobs['ext']) for f in item]
                dataDict[k]=Tframe

        # Version to loop over T and use setADMs
        # Better? More loops on IO, but simpler output!
        # CURRENTLY ONLY SUPPORTS pd.read_csv.
        if keyType == 'setADMsT':
            dataDict = {}
            for k,item in self.fileDict.items():

                Tframe = pd.concat([pd.read_csv(f, names=[f.split('_')[-1].strip(self.jobs['ext'])], header=None) for f in item], axis=1)
                dataDict[k]=Tframe

                Tcols = [f.split('_')[-1].strip(self.jobs['ext']) for f in item]

            dataDictT = {}
            for colName in Tcols:
                ADMs = []
                ADMLabels = []
                for k,df in dataDict.items():
                    ADMs.append(df[colName])   #*ADMscaleFactor)
                    ADMLabels.append([k[1],k[2],k[3]])

                # Add population (K=0) term?
                # This requires self.norm to be set at init.
                if addPop and hasattr(self,'norm'):
                    ADMs.append(np.ones(df[colName].size) * self.norm['K0'])
                    ADMLabels.append([0,0,0])

                dataDictT[colName] = setADMs(ADMs = ADMs,
                                                KQSLabels = np.array(ADMLabels),
                                                t=tIndex)  #, t=t)

                # TODO: add KQS to `checkSphDims` function?
                # See https://github.com/phockett/ePSproc/blob/a34745c2fdb8a2accc5a5fc71bbe14339ccf5e1c/epsproc/sphFuncs/sphConv.py#L194
#                 dataDictT[colName].attrs['harmonics'] = {'stackDim':{'ADM':['K','Q','S']}}

                # Additional renormalisation by (2*K+1)/8*pi^2
                # Generally shouldn't need this
                if renorm:
                    dataDictT[colName] = dataDictT[colName] * (2*dataDictT[colName].K+1)/(8*np.pi**2)

            self.dataDictT = dataDictT

            # Set also to standard self.data[key][dataType] style
            # This is used by existing core functionality
            for k in self.dataDictT.keys():
                self.data[k] = {'ADM':self.dataDictT[k]}

            # Stack to single dataarray
            # Use existing function, but note need to set coords curently
            Tkeys = list(self.data.keys())
            Tcoords = [int(item.rstrip('K')) for item in Tkeys]
            xrDA, xrDS, dataDict = datasetStack(self.data, dataType='ADM',
                                                stackDim = 'Temp', keys = Tkeys)
#             xrDA.assign_coords({"Temp":Tcoords})
            self.xr = xrDA.assign_coords({"Temp":Tcoords})
            self.xr.Temp.attrs['units'] = 'K'
            self.xr = self.xr.sortby('Temp')

        self.dataDict = dataDict
        self.TempIndex = Tcols # Assume all identical! BUT ORDERING MAY CHANGE!!!

        if self.verbose:
            print(f"Files from self.fileDict read OK, outputs: self.dataDict and self.data.")

    def setADMs(self,key=None,**kwargs):
        """
        Thin wrapper for ep.setADMs.

        Run `help(ep.setADMs)` for details.

        Additionally pass key=<str> to define data output, defaults to self.data['ADM'].

        """

        ADMs = setADMs(**kwargs)

        if key is None:
            key = 'ADM'

        self.data[key] = {'ADM':ADMs}

        if self.verbose:
            print(f"Set `self.data['{key}']['ADM']` from inputs.")


    def plot(self, keys = None, **kwargs):
        """
        Basic ADM plot with HV.
        See also ePSbase.ADMplot() for BLMplot() wrapper, although needs work.

        Parameters
        ----------
        keys : str or list, default = None
            Keys to use for plots.
            If None, use all keys.

        **kwargs : optional
            Additional args passes to hv.opts(**kwargs) for plot display control.

        TODO: tidy and return plot
        TODO: plot options
        TODO: dim stack options, currently just assumes Temp dim.
        TODO: selectors etc. Currently set for real data, all K,Q,S.
            (For single plot case K>0 set)

        """

        # Check keys, will set to all if keys=None
        keys = self._keysCheck(keys)

        # For single keys only...
        if len(keys) == 1:
#             ADMplot = self.dataDictT[keys[0]]  # Use dict?
            ADMplot = self.data[keys[0]]['ADM']  # Use data dict
            hvObj = ADMplot.unstack().where(ADMplot.unstack().K>0) \
                    .real.hvplot.line(x='t').overlay(['K','Q','S']).opts(width=700)
            showPlot(hvObj.opts(**kwargs), __notebook__=True)

        else:
        # ADMin.xr.unstack().squeeze().real.hvplot.line(x='t').overlay(['K'])
            keysT = [int(item.rstrip('K')) for item in keys]  # Convert to key coords

            hvObj = self.xr.sel(Temp=keysT).unstack().real.hvplot.line(x='t').overlay(['K','Q','S'])
            showPlot(hvObj.opts(**kwargs), __notebook__=True)

        #

    def normDict(self, normType = None):
        """
        List or set normalisation conventions (K_000 value).

        """

        norms = {}

        norms['sph'] = {'name':'Spherical Harmonics',
                        'K0': 1}
        norms['wignerDlinear'] = {'name':'Wigner D linear molecule',
                                  'K0': 1/(4*np.pi)}
        norms['wignerDpoly'] = {'name':'Wigner D polyatomic molecule',
                                'K0': 1/(8*np.pi**2)}

        self.norms = norms

        if normType is not None:
            self.norm = norms[normType]
            print(f"Set self.norm from self.norms['{normType}'].")


    def subsetADMs(self, dataKey = None, dataType = 'ADM', trange = None, tStep = 4):
        """
        Subselect ADMs to use for calcs.

        Basic case from https://phockett.github.io/ePSdata/OCS-preliminary/OCS_orbs8-11_AFBLMs_VM-ADMs_140122-JAKE_tidy-replot-200722_v5.html

        See also PEMtk fitting code for setSubset() method.

        """

        if dataKey is None:
            # Use first key for default case
            dataKey = self.TempIndex[0]

        if trange is None:
            # Set full axis, just downsample
            taxis = self.data[dataKey][dataType].t
            trange = [taxis[0], taxis[-1]]

        # Set ADMs to use
        ADMs = self.data[dataKey][dataType]

        # Selection & downsampling - adapted from https://epsproc.readthedocs.io/en/latest/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html#Test-compared-to-experimental-N2-AF-results...
        # See PEMtk for updated routines, https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full.html#Subselect-data
        # See Xarray docs for basics https://xarray.pydata.org/en/stable/user-guide/indexing.html#indexing-with-dimension-names

#         trange=[38, 44]  # Set range in ps for calc
#         tStep=2  # Set tStep for downsampling

        # SLICE version - was working, but not working July 2022, not sure if it's data types or Xarray version issue? Just get KeyErrors on slice.
        # 15/05/24: reimplemented this version... seems to be working, and avoids issues with dim changes in mask case.
        #           UPDATE: issues with older version of Xarray?
        #           Slice OK in xr2022, but fails in xr15.
        try:
            self.data['ADM'] = {'ADM': ADMs.sel(t=slice(trange[0],trange[1], tStep))}   # Set and update
        except KeyError:
            # Inds/mask version - seems more robust?
            # NOTE THIS ASSUMES DIMS!!!!
            tMask = (ADMs.t>trange[0]) & (ADMs.t<trange[1])
            # ind = np.nonzero(tMask)  #[0::tStep]
            # At = ADMs['time'][:,ind].squeeze()
            # ADMin = ADMs['ADM'][:,ind]

            self.data['ADM'] = {'ADM': ADMs[:,tMask][:,::tStep]}   # Set and update


        # TODO 15/05/24: need dim check and subselection here for dim change case too?
        if self.data['ADM']['ADM'].ndim > 2:
            print(f"*** Warning: setting ADMs with ndim = {self.data['ADM']['ADM'].ndim} may cause issues. Trying squeeze to fix...")
            self.data['ADM']['ADM'] = self.data['ADM']['ADM'].squeeze()

            if self.data['ADM']['ADM'].ndim < 3:
                print("Squeezed OK")
            else:
                print("Squeeze failed, additional subselection may be required.")


        # # Inds/mask version - seems more robust?
        # tMask = (ADMs.t>trange[0]) & (ADMs.t<trange[1])
        # # ind = np.nonzero(tMask)  #[0::tStep]
        # # At = ADMs['time'][:,ind].squeeze()
        # # ADMin = ADMs['ADM'][:,ind]
        #
        # self.data['ADM'] = {'ADM': ADMs[:,tMask][:,::tStep]}   # Set and update

        # Set metadata...
        self.data['ADM']['ADM'].attrs['subselection'] = {'trange':trange,
                                                     'tstep':tStep,
                                                     'sourceKey':dataKey,
                                                     'sourceDataType':dataType}

        print(f"Selecting {self.data['ADM']['ADM'].t.size} points")
        print("Set subset data to `self.data['ADM']['ADM']`")


        # NOTE: may want to apply more general methods per main plotting routines
        # E.G. from epsproc\classes\_plotters.py, line 772:
        # # 06/03/24 - reinstated Esubset for Erange settings.
        # # NEEDS TESTING - something here messes up padPlot() later for PL facetDims!!!
        # subset = self.Esubset(key = key, dataType = dataType, Etype = Etype, Erange = Erange)
        # subset = matEleSelector(subset, thres=thres, inds = selDims, dims = contiguousDims, sq = sqSelector)
