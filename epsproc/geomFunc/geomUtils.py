"""
ePSproc utility functions for geometric terms and tensors.

26/02/20    v1

Code from ePSproc_MFBLM_Numba_dev_tests_120220.py

"""


import numpy as np

# Basic function to match QN sets to NP array
# Input QNmask must be same dims as QNs, if set to None field will be skipped (i.e. match all values)
# Direct version
# %timeit QNsMask, mask = selQNsRow(QNs,[0, 0, 0], fields = [3,4,5], verbose = False)
# 19.7 µs ± 110 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
# With additional sorting, timings are similar
# %timeit QNsMask, mask = selQNsRow(QNs,[None, None, None, 0, 0, 0], verbose = False)
# 19.8 µs ± 59.5 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
def selQNsRow(QNs, QNmask, fields = None, verbose = 1):
    """
    Basic routine for selecting/filtering valid rows from an input array of QNs.

    This is similar to methods used in Matlab version, but likely rather slow.

    Parameters
    ----------
    QNs : np.array
        Values to be filtered.

    QNmask : list or np.array
        Values to be matched.
        Must be same dims as QNs, unless fields is also set.
        If set to None field will be skipped (i.e. match all values)

    fields : list or np.array, optional, default = None
        Which fields to match in QNs.

    verbose : int, optional, default = 1
        Print info to terminal.

    Returns
    -------
    QNs[boolMask] : np.array of selected values.

    boolMask : np.array corresponding to selected values.

    Examples
    ---------
    >>> QNsMask, mask = selQNsRow(QNs,[None, None, None, 0, 0, 0], verbose = False)
    >>> # With fields
    >>> QNsMask, mask = selQNsRow(QNs,[0, 0, 0], fields = [3,4,5], verbose = False)

    """

    # With sorting
    if fields is None:
        # # Sort input mask, reduce to masked fields only
        fields = []
        mask = []
        for n, item in enumerate(QNmask):
            if item is not None:
                fields.append(n)
                mask.append(item)

        fields = np.array(fields)
        mask = np.array(mask)
        boolMask = (QNs[:,fields] == mask).all(axis = 1)

    else:
        # Direct case with fields passed to function - timings are pretty much identical for test QNs, 423x6
        boolMask = (QNs[:,fields] == QNmask).all(axis = 1)

    if verbose == 1:
        print(f"{boolMask.sum()} rows selected, from {QNs.shape[0]} total {np.round(boolMask.sum()/QNs.shape[0]*100, 2)}%)")

    return QNs[boolMask], boolMask


# Generate 6D coords for Wigner 3j terms.
def genllL(Lmin = 0, Lmax = 10, mFlag = True):
    """
    Generate quantum numbers for angular momentum contractions (l, lp, L)

    Parameters
    ----------
    Lmin, Lmax : int, optional, default 0, 10
        Integer values for Lmin and Lmax respectively.

    mFlag : bool, optional, default = True
        m, mp take all values -l...+l if mFlag=True, or =0 only if mFlag=False

    Returns
    -------
    QNs : 2D np.array
        Values take all allowed combinations ['l','lp','L','m','mp','M'] up to l=lp=Lmax, one set per row.


    Examples
    ---------
    >>> # Calculate up to Lmax = 2
    >>> QNs = genllL(Lmax=2)
    >>> # Use with w3jTable function to calculate Wigner 3j terms
    >>> w3j = w3jTable(QNs = QNs)

    To do
    -----
    - Implement output options (see dev. function w3jTable).
    -
    """

    # Set QNs for calculation, (l,m,mp)
    QNs = []
    for l in np.arange(Lmin, Lmax+1):
        for lp in np.arange(Lmin, Lmax+1):

            if mFlag:
                mMax = l
                mpMax = lp
            else:
                mMax = 0
                mpMax = 0

            for m in np.arange(-mMax, mMax+1):
                for mp in np.arange(-mpMax, mpMax+1):
                    for L in np.arange(np.abs(l-lp), l+lp+1):
                        # Set M - note this implies specific phase choice.
                        # M = -(m+mp)
                        # M = (-m+mp)
                        # if np.abs(M) <= L:  # Skip terms with invalid M
                        #     QNs.append([l, lp, L, m, mp, M])

                        # Run for all possible M
                        for M in np.arange(-L, L+1):
                            QNs.append([l, lp, L, m, mp, M])

    return np.array(QNs)


# Generate 6D coords for Wigner 3j terms from a list of (l,l,L) coords.
def genllLList(Llist, uniqueFlag = True, mFlag = True):
    """
    Generate quantum numbers for angular momentum contractions (l, lp, L) from a passed list, (m, mp, M)=0 or all allowed terms.

    Parameters
    ----------
    Llist : list
        Values [l, lp, L] to use for calculations.
        Note this needs to be 2D array in current form of function, i.e. defined as np.array([[L1,L2,L3],...])

    uniqueFlag : bool, optional, default = True
        Drop duplicate [l,lp,L] sets from list.

    mFlag : bool, optional, default = True
        m, mp take all values -l...+l if mFlag=True, or =0 only if mFlag=False

    Returns
    -------
    QNs : 2D np.array
        Values take all allowed combinations ['l','lp','L','m','mp','M'] up to l=lp=Lmax, one set per row.


    Examples
    ---------
    >>> # Set from an array
    >>> QNs = genllLList(np.array([[1,1,2],[1,3,2],[1,1,2]]), mFlag = True)
    >>> # Use with w3jTable function to calculate Wigner 3j terms
    >>> w3j = w3jTable(QNs = QNs)

    To do
    -----
    - Implement output options (see dev. function w3jTable).
    -
    """

    # Select unique l sets
    if uniqueFlag:
        Llist = np.unique(Llist[:, 0:3], axis = 0)

    # Set QNs for calculation, (l,m,mp)
    # Slightly ugly repurposing of loop from genllL() code.
    QNs = []
    for lSet in Llist:
        l = lSet[0]
        lp= lSet[1]
        L = lSet[2]

        if mFlag:
            mMax = l
            mpMax = lp
        else:
            mMax = 0
            mpMax = 0

        for m in np.arange(-mMax, mMax+1):
            for mp in np.arange(-mpMax, mpMax+1):
                # for L in np.arange(np.abs(l-lp), l+lp+1):  # Removed for now, but may want to reinstate as valid L test?
                # Set M - note this implies specific phase choice.
                # M = -(m+mp)
                # M = (-m+mp)
                # if np.abs(M) <= L:  # Skip terms with invalid M
                #     QNs.append([l, lp, L, m, mp, M])

                # Run for all possible M
                for M in np.arange(-L, L+1):
                    QNs.append([l, lp, L, m, mp, M])

    return np.array(QNs)


# Generate 3j QNs lists from matrix elements (Xarray)
def genllpMatE(matE, uniqueFlag = True, mFlag = True, phaseConvention = 'S'):
    """
    Generate quantum numbers for angular momentum contractions (l, lp, L, m, mp, M) from sets of matrix elements.

    Parameters
    ----------
    matE : Xarray
        Xarray containing matrix elements, with QNs (l,m), as created by :py:func:`readMatEle`

    uniqueFlag : bool, default = True
        Check for duplicates and remove (can occur with some forms of matrix elements).

    mFlag : bool, optional, default = True
        m, mp take all passed values if mFlag=True, or =0 only if mFlag=False

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.
        If matE.attrs['phaseCons'] is already set, this will be used instead of passed args.


    Returns
    -------
    QNs : 2D np.array
        Values take all allowed combinations ['l','lp','L','m','mp','M'] from supplied matE


    Examples
    ---------
    >>> # From demo matrix elements
    >>> dataFile = os.path.join(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing
    >>> # Scan data file
    >>> dataSet = ep.readMatEle(fileIn = dataFile)
    >>> QNs = genllpMatE(dataSet[0])
    >>> # Use with w3jTable function to calculate Wigner 3j terms
    >>> w3j = w3jTable(QNs = QNs)

    To do
    -----
    - Implement output options (see dev. function w3jTable).
    -
    """

    # Local import.
    from epsproc.geomFunc.geomCalc import setPhaseConventions

    # For transparency/consistency with subfunctions, str/dict now set in setPhaseConventions()
    if 'phaseCons' in matE.attrs.keys():
        phaseCons = matE.attrs['phaseCons']
    else:
        phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # Get QNs from matE
    lList = matE.unstack().l.values  # Use unstack here, or np.unique(matE.l), to avoid duplicates

    # Use passed (m,mp) values, or run for m=mp=0 only.
    if mFlag:
        mList = matE.unstack().m.values
    else:
        mList = 0

    # Set QNs for calculation, (l,m,mp)
    QNs = []
    for l in lList:
        for lp in lList:
            for m in mList:
                for mp in mList:
                    for L in np.arange(np.abs(l-lp), l+lp+1):
                        # Set M - note this implies specific phase choice.
                        if phaseCons['genMatEcons']['negm']:
                            M = (-m+mp)  # Case for M -> -M switch
                        else:
                            M = -(m+mp)  # Usual phase convention.

                        # This is likely redundant/misguided, since already implied in phase convention above.
                        # if phaseCons['genMatEcons']['negM']:
                        #     M *= -1

                        if np.abs(M) <= L:  # Skip terms with invalid M
                            QNs.append([l, lp, L, m, mp, M])

                        # Run for all possible M
                        # for M in np.arange(-L, L+1):
                            # QNs.append([l, lp, L, m, mp, M])

    if uniqueFlag:
        return np.unique(QNs, axis = 0)
    else:
        return np.array(QNs)


#*************
# AF CASE Generate QNs - code adapted from MF cases above.
# - genKQSterms(): set as list of all possible values
# - genKQStermsFromTensors(): set from existing tensors

# Generate QNs for deltaKQS term - 3j product term
def genKQSterms(Kmin = 0, Kmax = 2, mFlag = True):
    # Set QNs for calculation, (l,m,mp)
    QNs = []
    for P in np.arange(0, 2+1):    # HARD-CODED R for testing - should get from EPR tensor defn. in full calcs.
        for K in np.arange(Kmin, Kmax+1):  # Eventually this will come from alignment term
            for L in np.arange(np.abs(P-K), P+K+1):  # Allowed L given P and K defined

                if mFlag:    # Include "m" (or equivalent) terms?
                    mMax = L
                    RMax = P
                    QMax = K
                else:
                    mMax = 0
                    RMax = 0
                    QMax = 0

                for R in np.arange(-RMax, RMax+1):
                    for Q in np.arange(-QMax, QMax+1):
                        #for M in np.arange(np.abs(l-lp), l+lp+1):
#                         for M in np.arange(-mMax, mMax+1):
                            # Set M - note this implies specific phase choice.
                            # M = -(m+mp)
                            # M = (-m+mp)
                            # if np.abs(M) <= L:  # Skip terms with invalid M
                            #     QNs.append([l, lp, L, m, mp, M])

                        # Run for all possible M
                        for M in np.arange(-L, L+1):
                            QNs.append([P, K, L, R, Q, M])

    return np.array(QNs)


# Generate QNs from EPR + AKQS tensors
def genKQStermsFromTensors(EPR, AKQS, uniqueFlag = True, phaseConvention = 'S'):
    '''
    Generate all QNs for :math:`\Delta_{L,M}(K,Q,S)` from existing tensors (Xarrays) :math:`E_{P,R}` and :math:`A^K_{Q,S}`.

    Cf. :py:func:`epsproc.geomFunc.genllpMatE`, code adapted from there.

    Parameters
    ----------
    matE : Xarray
        Xarray containing matrix elements, with QNs (l,m), as created by :py:func:`readMatEle`

    uniqueFlag : bool, default = True
        Check for duplicates and remove (can occur with some forms of matrix elements).

    mFlag : bool, optional, default = True
        m, mp take all passed values if mFlag=True, or =0 only if mFlag=False

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.
        If matE.attrs['phaseCons'] is already set, this will be used instead of passed args.


    Returns
    -------
    QNs1, QNs2 : two 2D np.arrays
        Values take all allowed combinations ['P','K','L','R','Q','M'] and ['P','K','L','Rp','S','S-Rp'] from supplied matE.
        Note phase conventions not applied to QN lists as yet.

    To do
    -----
    - Implement output options (see dev. function w3jTable).

    '''

    # Local import.
    from epsproc.geomFunc.geomCalc import setPhaseConventions

    # For transparency/consistency with subfunctions, str/dict now set in setPhaseConventions()
    if 'phaseCons' in EPR.attrs.keys():
        phaseCons = EPR.attrs['phaseCons']
    else:
        phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # Get QNs from inputs
    KQScoords = AKQS.unstack().coords  # Use unstack here, or np.unique(matE.l), to avoid duplicates
    PRcoords = EPR.unstack().coords

    # Use passed (m,mp) values, or run for m=mp=0 only.
#     if mFlag:
#         mList = matE.unstack().m.values
#     else:
#         mList = 0

    # Set QNs for calculation, one set for each 3j term
    QNs1 = []
    QNs2 = []
    for P in PRcoords['P'].values:   # Note dictionary syntax for coords objects
        for K in KQScoords['K'].values:
            for L in np.arange(np.abs(P-K), P+K+1):  # Allowed L given P and K defined

#                 if mFlag:    # Include "m" (or equivalent) terms?
#                     mMax = L
#                     RMax = P
#                     QMax = K
#                 else:
#                     mMax = 0
#                     RMax = 0
#                     QMax = 0

                for R in PRcoords['R'].values:
                    for Q in KQScoords['Q'].values:

                        # Set M, with +/- phase convention - TBC MAY BE INCORRECT IN THIS CASE/CONTEXT?
                        # Note that setting phaseCons['afblmCons']['negM']  = phaseCons['genMatEcons']['negm'] is current default case, but doesn't have to be!
                        if phaseCons['genMatEcons']['negm']:
                            M = (-R+Q)  # Case for M -> -M switch
                        else:
                            M = -(R+Q)  # Usual phase convention.

                        if (abs(R)<=P) and (abs(Q)<=K) and (abs(M)<=L): # Check term is valid.
                            QNs1.append([P, K, L, R, Q, M])

                        # QNs1.append([P, K, L, R, Q, M])

                # Set Rp and S - these are essentially independent of R,Q,M, but keep nested for full dim output.
                # for Rp in PRcoords['P'].values:
                # for Rp in PRcoords['R'].values:  # OOOPS, this gives accidental correlations between R and Rp
                for Rp in np.arange(-P, P+1):
                    for S in KQScoords['S'].values:
                        # SRp = S-Rp  # Set final 3j term, S-Rp
                        if phaseCons['genMatEcons']['negm']:
                            SRp = (-Rp+S)  # Case for M -> -M switch
                        else:
                            SRp = -(Rp+S)  # Usual phase convention.
#                             SRp = S-Rp  # Set final 3j term, S-Rp
                            # SRp = -(Rp-S)  # Usual phase convention.

                        if (abs(Rp)<=P) and (abs(S)<=K) and (abs(SRp)<=L): # Check term is valid.
                            QNs2.append([P, K, L, Rp, S, SRp])

                        # QNs2.append([P, K, L, Rp, S, SRp])

                            #for M in np.arange(np.abs(l-lp), l+lp+1):
    #                         for M in np.arange(-mMax, mMax+1):
                                # Set M - note this implies specific phase choice.
                                # M = -(m+mp)
                                # M = (-m+mp)
                                # if np.abs(M) <= L:  # Skip terms with invalid M
                                #     QNs.append([l, lp, L, m, mp, M])

                            # Run for all possible M
    #                         for M in np.arange(-L, L+1):
    #                             QNs.append([P, K, L, R, Q, M])


    if uniqueFlag:
        return np.unique(QNs1, axis = 0), np.unique(QNs2, axis = 0)
    else:
        return np.array(QNs1), np.array(QNs2)
