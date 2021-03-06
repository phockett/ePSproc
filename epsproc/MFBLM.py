r"""
ePSproc MFBLM functions
-----------------------

23/09/19        Testing function caching...

19/09/19        Working & verified basic version (slow loop).

05/09/19    v1  Initial python version.
                Based on original Matlab code ePS_MFBLM.m


Formalism
---------

Should be something like this, with possible substitutions or phase swaps.

.. math::
    \begin{eqnarray}
    \beta_{L,-M}^{\mu_{i},\mu_{f}} & = & \sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{M}(-1)^{m}(-1)^{(\mu'-\mu_{0})}\left(\frac{(2l+1)(2l'+1)(2L+1)}{4\pi}\right)^{1/2}\left(\begin{array}{ccc}
    l & l' & L\\
    0 & 0 & 0
    \end{array}\right)\left(\begin{array}{ccc}
    l & l' & L\\
    -m & m' & -M
    \end{array}\right)\nonumber \\
     & \times & \sum_{P,R',R}(2P+1)(-1)^{(R'-R)}\left(\begin{array}{ccc}
    1 & 1 & P\\
    \mu & -\mu' & R'
    \end{array}\right)\left(\begin{array}{ccc}
    1 & 1 & P\\
    \mu_{0} & -\mu_{0} & R
    \end{array}\right)D_{-R',-R}^{P}(R_{\hat{n}})I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)
    \end{eqnarray}


Exact numerics may vary.


"""

# Imports
import numpy as np
import pandas as pd
import xarray as xr
from functools import lru_cache  # For function result caching

# Special functions
# from scipy.special import sph_harm
import spherical_functions as sf
import quaternion

# Package fns.
# TODO: tidy this up!
from epsproc.util import matEleSelector


# Wrappers for function caching
# Tested for mfblm calcs, N2 vs. E, NO2 vs. Euler angs. (in test Jupyter notebook).  No change vs. bare functions in this case.
@lru_cache(maxsize = None)
def Wigner3jCached(j_1,j_2,j_3,m_1,m_2,m_3):
    """Wrapper for 3j caching with functools.lru_cache"""
    return sf.Wigner3j(j_1,j_2,j_3,m_1,m_2,m_3)

@lru_cache(maxsize = None)
def Wigner_D_element_Cached(pRot, P, Rp, R):
    """Wrapper for WignerD caching with functools.lru_cache"""
    return  sf.Wigner_D_element(pRot, P, Rp, R)

# Create Xarray from set of BLMs
def blmXarray(BLM, Eke):
    """
    Create Xarray from BLM list, format BLM = [L, M, Beta], at a single energy Eke.

    Array sorting only valid for 2D BLM array, for B00=0 case pass BLM = None.

    """

    # Set up Xarray
    if BLM is not None:
        QNs = pd.MultiIndex.from_arrays(BLM[:,0:2].real.T.astype('int8'), names = ['l','m'])  # Set small (l,m) here for compatibility with sphCalc()
        BLMX = xr.DataArray(BLM[:,2], coords={'BLM':QNs}, dims = ['BLM'])
        # Set other dims - this might fail for multi symmetry cases... or if pre-selected?  Best done in calling/chunking function.
        # BLMX = BLMX.expand_dims({'Sym':[matE.Sym.data], 'Eke':[matE.Eke.data]})

    # Case for BLM = None, set array with B00 = 0
    else:
        QNs = pd.MultiIndex.from_arrays([[0],[0]], names = ['l','m'])  # Set small (l,m) here for compatibility with sphCalc()
        BLMX = xr.DataArray([0], coords={'BLM':QNs}, dims = ['BLM'])


    BLMX = BLMX.expand_dims({'Eke':[Eke]})

    BLMX.attrs['dataType'] = 'BLM'

    return BLMX


# Try basic looping version, similar to Matlab code...
def MFBLMCalcLoop(matE, eAngs = [0,0,0], thres = 1e-6, p=0, R=0, verbose=1):
    """
    Calculate inner loop for MFBLMs, based on passed set of matrix elements (Xarray).

    Loop code based on Matlab code ePSproc_MFBLM.m

    This works, but not very clean or efficient - should be possible to parallelize & vectorise... but this is OK for testing, benchmarking and verification purposes.

    Parameters
    ----------
    matE : Xarray
        Contains one set of matrix elements to use for calculation.
        Currently assumes these are a 1D list, with associated (l,m,mu) parameters, as set by :py:func:`epsproc.MFBLM.mfblm()`.

    eAngs : [phi,theta,chi], optional, default = [0,0,0]
        Single set of Euler angles defining polarization geometry.

    thres : float, optional, default = 1e-4
        Threshold value for testing significance of terms. Terms < thres will be dropped.

    p : int, optional, default p = 0
        LF polarization term. Currently only valid for p = 0

    R : int, optional, default R = 0
        LF polarization term (from tensor contraction). Currently only valid for p = 0

    verbose : int, optional, default 1
        Verbosity level:

         - 0: Silent run.
         - 1: Print basic info.
         - 2: Print intermediate C parameter array to terminal when running.

    Returns
    -------
    BLMX : Xarray
        Set of B(L,M; eAngs, Eke) terms for supplied matrix elements, in an Xarray.
        For cases where no values are calculated (below threshold), return an array with B00 = 0 only.

    Limitations \& To Do
    --------------------

    * Currently set with (p,R) values passed, but only valid for (0,0) (not full sum over R terms as shown in formalism above.)
    * Set to accept a single set of matrix elements (single E), assuming looping over E (and other parameters) elsewhere.
    * Not explicitly parallelized here, should be done by calling function.
    (Either via Xarray methods, or numba/dask...?  http://xarray.pydata.org/en/stable/computation.html#wrapping-custom-computation)
    * Coded for ease, not efficiency - there will be lots of repeated ang. mom. calcs. when run over many sets of matrix elements.
    .. rst-class:: strike
    * Scale factor currently not propagated.

    * (SF now propagated via Xarrays and implemented in main calling function )

    """

    # Issues with this, since I can't seem to get coords out once at singleton dim level (issue with pd multilevel?)
#    for matE1 in matE:
#        for matE2 in matE:
#            for L in np.arange(0, matE1.LM[0] + matE2.l +1):
#                print(L)


    # Generate LM pairs list via indexing
    LMlist = pd.MultiIndex.from_product([matE.SumDim, matE.SumDim], names = ['LM1','LM2'])
    indList = pd.MultiIndex.from_product([np.arange(0, matE.size), np.arange(0, matE.size)], names = ['ind1','ind2'])

    if verbose > 0:
        print('Calculating MFBLMs for {0} pairs... Eke = {1} eV, eAngs = ({2})'.format(indList.size, matE.Eke.data, eAngs))

    # Set pol terms - now passed
    # p = 0
    # R = 0 # For single pol state only, i.e. p-p'=0
    # pRot = quaternion.from_euler_angles(alpha, beta, gamma)
    # pRot = quaternion.from_euler_angles(0,0,0)
    if not isinstance(eAngs, xr.DataArray):
        pRot = quaternion.from_euler_angles(eAngs)
    else:
        pRot = eAngs.values.item()  # For Xarray case will already be passed as quaternion. Get issues later if value in np array.

    # Set variable for outputs
    C = []

    # Loop over all pairs defined in lists & calculate contributions
    for rows in indList:
       # print(matE[n[0]], matE[n[1]])

       # For clarity, set values explicitly here
       # Explicit casting to int() may also be necessary in some cases (with passing of Xarray values to WignerD/3j sometimes get type issues)
       # In usual notation, r1 == unprimed terms, r2 == primed terms
       r1 = rows[0]
       r2 = rows[1]

       l1 = int(matE.l[r1])
       m1 = int(matE.m[r1])
       # mu1 = int(matE[r1].mu)  # Case for single mu
       mu1 = int(matE.mu[r1])  # Case for stacked dataArray (mu per item)

       l2 = int(matE.l[r2])
       m2 = int(matE.m[r2])
       # mu2 = int(matE[r2].mu) # Case for single mu
       mu2 = int(matE.mu[r2])  # Case for stacked dataArray (mu per item)

       matEprod = matE[r1] * matE[r2].conj()   # Mat element product
                                               #TODO: could set additional thresholding here for speed

       # Set allowed M value, note phase convention (tested in Matlab version)
       M = (-m1 + m2)

       if np.abs(matEprod) > thres:

           # Gamma terms for polarization - sum over P
           gammaP = 0

           # Loop over allowed P = [0...2] for 2-photon case
           for P in np.arange(0, 3):
               Rp = mu2 - mu1     # mu2-mu1, includes Rp > -Rp for numerics.

               if np.abs(Rp) <= P:  # Check for allowed terms

                   # Sum over R,R' projections omitted - only R=0 set above
                   # Note polarization terms, here mu is MF, and p is LF.
                   # Note use of matEprod.data here - otherwise output is kept as Xarray (which could be useful in future...)
                   gammaP = gammaP + (2*P+1) * (-1)**(Rp-R) * Wigner3jCached(1,1,P,mu1,-mu2,Rp) * Wigner3jCached(1,1,P,p,-p,R) * Wigner_D_element_Cached(pRot, P,-Rp,-R).conj()

           # Loop over allowed B(LM) and calculate
           # NOTE - for multiindex coords matE.l[r1] works, but matE[r1].l DOESN'T (Xarray v0.12.3)
           if np.abs(gammaP) > thres:
               for L in np.arange(0, l1 + l2 +1):
                   # print(L)

                   # Calculate associated gamma term, (L,M) values
                   gammaLM = Wigner3jCached(l1, l2, L, 0, 0, 0) * Wigner3jCached(l1, l2, L, -m1, m2, -M)

                   # Set any remaining terms and sum
                   phase = (-1)**M *(-1)**m1 * (-1)**(mu2-p)
                   degen = np.sqrt(((2*l1+1)*(2*l2+1)*(2*L+1))/(4*np.pi))

                   gamma = gammaLM*gammaP*phase*degen

                   LMterm = gamma*matEprod.data   # *(sf^2)  #TODO: propagate SF

                   C.append([l1, m1, l2, m2, L, M, LMterm, gammaLM, gammaP, phase, degen])
                   # print(C[-1])


    if len(C) > 0:
        # Threshold terms
        C = np.array(C)
        # print(C)
        indThres = np.abs(C[:,6]) > thres
        Cthres = C[indThres,:]

        # Sun values that remain after thresholding
        if Cthres.size > 0:
            if verbose > 1:
                print(Cthres)

            # Calculate output values as sum over (L,M) terms
            BLM = []
            for L in np.arange(0, int(Cthres[:,4].max().real)+1):
               for M in np.arange(-L, L+1):
                   maskLM = (Cthres[:,4].real.astype(int)==L)*(Cthres[:,5].real.astype(int)==M);
                   BLM.append([L, M, Cthres[maskLM,6].sum()])

            BLM = np.array(BLM)

            # # Set output as Xarray with coords - NOW moved to function
            # QNs = pd.MultiIndex.from_arrays(BLM[:,0:2].real.T.astype('int8'), names = ['l','m'])  # Set small (l,m) here for compatibility with sphCalc()
            # BLMX = xr.DataArray(BLM[:,2], coords={'BLM':QNs}, dims = ['BLM'])
            # # Set other dims - this might fail for multi symmetry cases... or if pre-selected?  Best done in calling/chunking function.
            # # BLMX = BLMX.expand_dims({'Sym':[matE.Sym.data], 'Eke':[matE.Eke.data]})
            # BLMX = BLMX.expand_dims({'Eke':[matE.Eke]})

            # Set output as Xarray with coords
            BLMX = blmXarray(BLM, matE.Eke)

        else:
            BLMX = blmXarray(None, matE.Eke)

    else:
        # BLMX = None
        BLMX = blmXarray(None, matE.Eke)

    return BLMX


# Master BLM calculation routine.
# TODO: set for differnt calc routines, once developed!
def mfblm(daIn, selDims = {'Type':'L'}, eAngs = [0,0,0], thres = 1e-4, sumDims = ('l','m','mu','Cont','Targ','Total','it'), SFflag = True, verbose = 1):
    """
    Calculate MFBLMs for a range of (E, sym) cases. Default is to calculated for all symmetries at each energy.

    Parameters
    ----------
    da : Xarray
        Contains matrix elements to use for calculation.
        Matrix elements will be sorted by energy and BLMs calculated for each set.

    selDims : dict, optional, default = {'Type':'L'}
        Additional sub-selection to be applied to matrix elements before BLM calculation.
        Default selects just length gauge results.

    eAngs : [phi,theta,chi], optional, default = [0,0,0]
        Single set of Euler angles defining polarization geometry.

    thres : float, optional, default = 1e-4
        Threshold value for testing significance of terms. Terms < thres will be dropped.

    sumDims : tuple, optional, default = ('l','m','mu','Cont','Targ','Total','it')
        Defines which terms are summed over (coherent) in the MFBLM calculation.
        (These are used to flatten the Xarray before calculation.)
        Default includes sum over (l,m), symmetries and degeneracies (but not energies).

    SFflag : bool, default = True
        Normalise by scale factor to give X-sections (B00) in Mb

    verbose : int, optional, default 1
        Verbosity level:

         - 0: Silent run.
         - 1: Print basic info.
         - 2: Print intermediate C parameter array to terminal when running.

    Returns
    -------
    Xarray
        Calculation results BLM, dims (Euler, Eke, l,m). Some global attributes are also appended.

    Limitations
    -----------
    Currently set to loop calcualtions over energy only, and all symmetries.
    Pass single {'Cont':'sym'} to calculated for only one symmetry group.

    TODO: In future this will be more elegant.
    TODO: Setting selDims in output structure needs more thought for netCDF save compatibility.

    """

    # Make explicit copy of data to avoid any overwrite issues
    da = daIn.copy()
    da.attrs = daIn.attrs

    # Use SF (scale factor)
    # Write to data.values to make sure attribs are maintained. (Not the case for da = da*da.SF)
    if SFflag:
        da.values = da * da.SF

    # Unstack & sub-select data array
    daUnStack = da.unstack()
    daUnStack = matEleSelector(daUnStack, thres = thres, inds = selDims, sq = True)

    # Check dims vs. inds selection and/or key dims in output
    for indTest in sumDims:
        # if (indTest in inds) or (indTest not in unstackTest.dims):
        if indTest not in daUnStack.dims:
            daUnStack = daUnStack.expand_dims(indTest)

    # Restack array along sumInds dimensions, drop NaNs and squeeze.
    daSumDim = daUnStack.stack(SumDim = sumDims).dropna('SumDim').squeeze()

    # Check dims, and loop over any other dims.
    # TODO: see http://xarray.pydata.org/en/stable/groupby.html#multidimensional-grouping
    # For now assume Eke loop only, since symmetries already stacked above - will need to remove that to be more general here.
    # dimList = []
    # # [dimList.append(indTest) if indTest is not 'SumDim' for indTest in daSumDim.dims]
    # for indTest in daSumDim.dims:
    #     if indTest is not 'SumDim':
    #         dimList.append(indTest)

    # Loop over Eke and calculate BLMs
    BLMXlist = []
    if daSumDim.Eke.size > 1:
        for n, daE in daSumDim.groupby('Eke'):
            if SFflag:
                BLMXlist.append(MFBLMCalcLoop(daE, eAngs = eAngs, thres = thres, verbose = verbose)) #/daE.SF)  # Rescale by SF if required
            else:
                BLMXlist.append(MFBLMCalcLoop(daE, eAngs = eAngs, thres = thres, verbose = verbose))

            # Add dims - currently set for Euler angles only (single set)
            # Can't seem to add mutiindex as a single element, so set dummy coord here and replace below.
            # BLMXlist[-1] = BLMXlist[-1].expand_dims({'Euler':[n]})  # Set as index

        # Restack along Eke
        BLMXout = xr.combine_nested(BLMXlist, concat_dim=['Eke'])

    else:
        BLMXout = MFBLMCalcLoop(daSumDim, eAngs = eAngs, thres = thres)  # .expand_dims({'Euler':[0]})

    # Set Euler angles used (cf. MFPAD.py code), assumed same for each Eke
    # May fail for singleton Eke dim...?
    # Euler = pd.MultiIndex.from_arrays(np.tile(eAngs,(BLMXout.Eke.size,1)).T, names = ['P','T','C'])
    # BLMXout = BLMXout.assign_coords(Euler = Euler)

    if isinstance(eAngs, xr.DataArray):
        # For Xarray defined case, use existing multindex & propagate corresponding labels.
        # BLMXout = BLMXout.expand_dims({'Euler':eAngs.get_index('Euler')})  # THIS DOESN'T WORK FOR A 0D ARRAY, DOH.
        Euler = pd.MultiIndex.from_arrays(np.tile(np.asarray(eAngs.Euler.values.item()),(1,1)).T, names = ['P','T','C'])    # Reconstruct as for normal case below.
        BLMXout = BLMXout.expand_dims({'Euler':Euler})
        # BLMXout['Labels']=('Euler', eAngs.Labels.values)  # Now set for all eAngs in mfblmEuler
    else:
        # For other cases set from eAngs array values directly.
        # Set instead as singleton dim
        Euler = pd.MultiIndex.from_arrays(np.tile(eAngs,(1,1)).T, names = ['P','T','C'])    # Left tile code here to prevent pd errors on lists.
        BLMXout = BLMXout.expand_dims({'Euler':Euler})


    # Fix XS issue due to SF^2 - in tests this matchs GetCro results
    # ISSUE: in current form, with SF(Eke,Sym), this reintroduces Sym axis even if already summed over.
    # Q: Is SF sym dependent...?  If not, remove link.  If so, sum before division.
    # NOW: checked in dumpIdySegsParseX, set to single dim (Eke only) if diff < machine epsilon
    if SFflag:
        # Check dims on SF are OK...
        if da.SF.ndim > 1:
            print(f'*** Warning: SF has {0} dims, skipping renormalisation of BLMs by 1/SF.', da.SF.ndim)
        else:
            BLMXout = BLMXout / da.SF  # Renorm - should be able to sort this so it's not required...?  Multiply MFBLMCalcLoop result by SF only, rather than all matE?

    # Reorder & normalise by X-sect (B00)
    BLMXout = BLMXout.transpose()
    BLMXout['XS'] = (('Eke','Euler'), BLMXout[0].data)  # Set XS = B00
    BLMXout = BLMXout/BLMXout.XS  # Normalise

    # Set/propagate global properties
    BLMXout.attrs = da.attrs
    BLMXout.attrs['thres'] = thres
    BLMXout.attrs['sumDims'] = sumDims # May want to explicitly propagate symmetries here...?
    BLMXout.attrs['selDims'] = [(k,v) for k,v in selDims.items()]  # Can't use Xarray to_netcdf with dict set here, at least for netCDF3 defaults.
    BLMXout.attrs['dataType'] = 'BLM'

    return BLMXout


# Wrapper for mfblm for a set of Euler angles
def mfblmEuler(da, selDims = {'Type':'L'}, eAngs = [0,0,0], thres = 1e-4, sumDims = ('l','m','mu','Cont','Targ','Total','it'), SFflag = True, verbose = 1):
    """
    Wrapper for :py:func:`epsproc.MFBLM.mfblm()` for a set of Euler angles.  All other parameters are simply passed to mfblm().
    Calculate MFBLMs for a range of (E, sym) cases. Default is to calculated for all symmetries at each energy.

    Parameters
    ----------
    da : Xarray
        Contains matrix elements to use for calculation.
        Matrix elements will be sorted by energy and BLMs calculated for each set.

    selDims : dict, optional, default = {'Type':'L'}
        Additional sub-selection to be applied to matrix elements before BLM calculation.
        Default selects just length gauge results.

    eAngs : [phi,theta,chi], optional, default = [0,0,0]
        Set of Euler angles defining polarization geometry.
        List or np.array, dims(N, 3).
        Basic Xarray support also in place, as generated by :py:func:`epsproc.sphCalc.setPolGeoms()`

    thres : float, optional, default = 1e-4
        Threshold value for testing significance of terms. Terms < thres will be dropped.

    sumDims : tuple, optional, default = ('l','m','mu','Cont','Targ','Total','it')
        Defines which terms are summed over (coherent) in the MFBLM calculation.
        (These are used to flatten the Xarray before calculation.)
        Default includes sum over (l,m), symmetries and degeneracies (but not energies).

    SFflag : bool, default = True
        Normalise by scale factor to give X-sections (B00) in Mb

    verbose : int, optional, default 1
        Verbosity level:

         - 0: Silent run.
         - 1: Print basic info.
         - 2: Print intermediate C parameter array to terminal when running.

    Returns
    -------
    Xarray
        Calculation results BLM, dims (Euler(P,T,C), Eke, BLM(l,m)) - as per :py:func:`epsproc.util.BLMdimList`. Some global attributes are also appended.

    Limitations
    -----------
    Currently set to loop calcualtions over energy only, and all symmetries.
    Pass single {'Cont':'sym'} to calculated for only one symmetry group.

    TODO: In future this will be more elegant.
    TODO: Setting selDims in output structure needs more thought for netCDF save compatibility.

    """

    # Check size of passed set of Euler angles
    # For ease of manipulation, just change to np.array if necessary!
    if isinstance(eAngs, list):
        eAngs = np.array(eAngs)

    if isinstance(eAngs, np.ndarray):
        # For a single set of eAngs, just pass directly.
        if eAngs.ndim == 1:
            BLMXout = mfblm(da, selDims = selDims, eAngs = eAngs, thres = thres, sumDims = sumDims, SFflag = SFflag, verbose = verbose)
        else:
        # Loop over eAngs and calculate
            BLM = []
            for angsIn in range(0, eAngs.shape[0]):
                BLM.append(mfblm(da, selDims = selDims, eAngs = eAngs[angsIn,:], thres = thres, SFflag = SFflag, verbose = verbose))

            # Stack results
            BLMXout = xr.combine_nested(BLM,'Euler')

    # Quick hack for Xarray case - 1D array with quaternions already defined, as returned by `setPolGeoms()`
    if isinstance(eAngs, xr.DataArray):
        BLM = []
        for angsIn in range(0, eAngs.shape[0]):
            BLM.append(mfblm(da, selDims = selDims, eAngs = eAngs[angsIn], thres = thres, SFflag = SFflag, verbose = verbose))

        # Stack results
        BLMXout = xr.combine_nested(BLM,'Euler')
        BLMXout['Labels']=('Euler', eAngs.Labels.values) # Propagate labels

    return BLMXout
