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
                        M = -(m+mp)
                        if np.abs(M) <= L:  # Skip terms with invalid M
                            QNs.append([l, lp, L, m, mp, M])

    return np.array(QNs)


# Generate 3j QNs lists from matrix elements
def genllpMatE(matE, uniqueFlag = True, mFlag = True):
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
                        M = -(m+mp)
                        if np.abs(M) <= L:  # Skip terms with invalid M
                            QNs.append([l, lp, L, m, mp, M])

    if uniqueFlag:
        return np.unique(QNs, axis = 0)
    else:
        return np.array(QNs)
