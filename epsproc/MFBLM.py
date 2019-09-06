r"""
ePSproc MFBLM functions
-----------------------


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
# Special functions
# from scipy.special import sph_harm
import spherical_functions as sf
import quaternion


# Try basic looping version, similar to Matlab code...
def MFBLMCalcLoop(matE, eAngs = [0,0,0], thres = 1e-10, p=0, R=0):
    """
    Calculate inner loop for MFBLMs, based on passed set of matrix elements (Xarray).

    Loop code based on Matlab code ePSproc_MFBLM.m

    This works, but not very clean or efficient - should be possible to parallelize & vectorise...

    Parameters
    ----------
    matE : Xarray
        Contains one set of matrix elements to use for calculation.

    eAngs : [phi,theta,chi], optional, default = [0,0,0]
        Single set of Euler angles defining polarization geometry.

    thres : float, optional, default = 1e-10
        Threshold value for testing significance of terms. Terms < thres will be dropped.

    p : int, optional, default p = 0
        LF polarization term. Currently only valid for p = 0

    R : int, optional, default R = 0
        LF polarization term (from tensor contraction). Currently only valid for p = 0

    Returns
    -------
    BLMX : Xarray
        Set of B(L,M; eAngs, Eke) terms for supplied matrix elements, in an Xarray.

    Limitations & To Do
    --------------------

    * Currently set with (p,R) values passed, but only valid for (0,0) (not full sum over R terms as shown in formalism above.)
    * Set to accept a single set of matrix elements (single E), assuming looping over E (and other parameters) elsewhere.
    * Not explicitly parallelized here, should be done by calling function.
    * Coded for ease, not efficiency - there will be lots of repeated ang. mom. calcs. when run over many sets of matrix elements.
    * Scale factor currently not propagated.

    """

    # Issues with this, since I can't seem to get coords out once at singleton dim level (issue with pd multilevel?)
#    for matE1 in matE:
#        for matE2 in matE:
#            for L in np.arange(0, matE1.LM[0] + matE2.l +1):
#                print(L)


    # Generate LM pairs list via indexing
    LMlist = pd.MultiIndex.from_product([matE.LM, matE.LM], names = ['LM1','LM2'])
    indList = pd.MultiIndex.from_product([np.arange(0, matE.size), np.arange(0, matE.size)], names = ['ind1','ind2'])

    print('Calculating MFBLMs for {0} pairs...'.format(indList.size))

    # Set pol terms - now passed
    # p = 0
    # R = 0 # For single pol state only, i.e. p-p'=0
    # pRot = quaternion.from_euler_angles(alpha, beta, gamma)
    # pRot = quaternion.from_euler_angles(0,0,0)
    pRot = quaternion.from_euler_angles(eAngs)

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
       mu1 = int(matE[r1].mu)
       l2 = int(matE.l[r2])
       m2 = int(matE.m[r2])
       mu2 = int(matE[r2].mu)

       # Loop over allowed B(LM) and calculate
       # NOTE - for multiindex coords matE.l[r1] works, but matE[r1].l DOESN'T (Xarray v0.12.3)
       for L in np.arange(0, l1 + l2 +1):
           # print(L)

           # Set allowed M value, note phase convention (tested in Matlab version)
           M = (-m1 + m2)

           # Calculate associated gamma term, (L,M) values
           gammaLM = sf.Wigner3j(l1, l2, L, 0, 0, 0) * sf.Wigner3j(l1, l2, L, -m1, m2, -M)

           # Gamma terms for polarization - sum over P
           gammaP = 0
           matEprod = matE[r1] * matE[r2].conj()   # Mat element product
                                                   #TODO: could set additional thresholding here for speed

            # Loop over allowed P = [0...2] for 2-photon case
           for P in np.arange(0, 3):
               Rp = mu2 - mu1     # mu2-mu1, includes Rp > -Rp for numerics.

               if np.abs(Rp) <= P:  # Check for allowed terms

                   # Sum over R,R' projections omitted - only R=0 set above
                   # Note polarization terms, here mu is MF, and p is LF.
                   # Note use of matEprod.data here - otherwise output is kept as Xarray (which could be useful in future...)
                   gammaP = gammaP + (2*P+1) * (-1)**(Rp-R) * sf.Wigner3j(1,1,P,mu1,-mu2,Rp) * sf.Wigner3j(1,1,P,p,-p,R) * sf.Wigner_D_element(pRot, P,-Rp,-R).conj() * matEprod.data


           # Set any remaining terms and sum
           phase = (-1)**M *(-1)**m1 * (-1)**(mu2-p)
           degen = np.sqrt(((2*l1+1)*(2*l2+1)*(2*L+1))/(4*np.pi))

           gamma = gammaLM*gammaP*phase*degen

           LMterm = gamma   # *(sf^2)  #TODO: propagate SF

           C.append([l1, m1, l2, m2, L, M, LMterm, gammaLM, gammaP, phase, degen])
           # print(C[-1])


    # Threshold terms
    C = np.array(C)
    indThres = np.abs(C[:,7]) > thres
    Cthres = C[indThres,:]

    # Calculate output values as sum over (L,M) terms
    BLM = []
    for L in np.arange(0, int(Cthres[:,5].max().real)+1):
       for M in np.arange(-L, L+1):
           maskLM = (Cthres[:,5].real.astype(int)==L)*(Cthres[:,6].real.astype(int)==M);
           BLM.append([L, M, Cthres[maskLM,7].sum()])

    BLM = np.array(BLM)

    # Set output as Xarray with coords
    QNs = pd.MultiIndex.from_arrays(BLM[:,0:2].real.T.astype('int8'), names = ['L','M'])
    BLMX = xr.DataArray(BLM[:,2], coords={'BLM':QNs}, dims = ['BLM'])
    # BLMX = BLMX.expand_dims({'Sym':[matE.Sym.data], 'Eke':[matE.Eke.data]})           # This might fail for multi symmetry cases... or if pre-selected?
    BLMX = BLMX.expand_dims({'Eke':[matE.Eke]})

    return BLMX
