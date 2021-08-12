"""
ePSproc Wigner 3,6,9j implementation with Sympy

Basic wrappers to tabulate Sympy results.
For preference, use w3jVecMethods which contain functions for accelerated (vector, parallel, GPU) Wigner 3j calculations with spherical_functions library.

07/11/20    v1  Adding for testing half-integer values, and 6j capabilities.

"""


try:
    import sympy as sp
    from sympy.physics.wigner import wigner_3j, wigner_6j, wigner_9j

except ImportError as e:
    if e.msg != "No module named 'sympy'":
        raise
    print('* Sympy not found, Wigner Nj functions not available.')


# Bare function with loop, pass w3j_QNs
def w3jSympy(QNs, Nflag = True):
    """
    Basic wrapper for Sympy Wigner 3j function for a list of QNs.

    For function details, see the Sympy documentation, https://docs.sympy.org/latest/modules/physics/wigner.html

    For a higher-level wrapper (recommended), use :py:func:`epsproc.geomFunc.geomCalc.w3jTable()`, with backend = 'sympy'.

    Parameters
    ----------

    QNs : list
        List of QNs [l, lp, L, m, mp, M] to compute 3j terms for.

    Nflag : bool, optional, default = True
        Return float or symbolic values.

    Returns
    -------

    w3j_QNs : list
        List of w3j values.

    """
    # w3j_QNs = np.zeros(QNs.shape[0])
    # for n in range(QNs.shape[0]):
    #     w3j_QNs[n] = sf.Wigner3j(QNs[n,0], QNs[n,1], QNs[n,2], QNs[n,3], QNs[n,4], QNs[n,5])

    if Nflag:
        # Convert to sypmy.core.float, then force to native python float 
        w3j_QNs = [float(sp.N(wigner_3j(QN[0], QN[1], QN[2], QN[3], QN[4], QN[5]))) for QN in QNs]
    else:
        w3j_QNs = [wigner_3j(QN[0], QN[1], QN[2], QN[3], QN[4], QN[5]) for QN in QNs]

    return w3j_QNs
