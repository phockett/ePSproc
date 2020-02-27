"""
ePSproc Wigner 3j vector methods

Functions for accelerated (vector, parallel, GPU) Wigner 3j calculations.

26/02/20    v1  Functions for accelerated Wigner 3j calcs plus sorting.

"""

# Imports
import numpy as np
import numba as nb
import spherical_functions as sf

#*****************  Wrappers for functions
# For tests, see
# Test wrap sf.Wigner3j in Numba

# Bare function with loop, pass w3j_QNs
def Wigner3jQNs(QNs):
    w3j_QNs = np.zeros(QNs.shape[0])
    for n in range(QNs.shape[0]):
        w3j_QNs[n] = sf.Wigner3j(QNs[n,0], QNs[n,1], QNs[n,2], QNs[n,3], QNs[n,4], QNs[n,5])
    return w3j_QNs

# Cached function, only for single set of QNs.
# @lru_cache(maxsize = None)
# def Wigner3jLRUCached(QNlist):
#     """Wrapper for 3j caching with functools.lru_cache"""
#     return sf.Wigner3j(QNlist[0], QNlist[1], QNlist[2], QNlist[3], QNlist[4], QNlist[5])

# Try joblib Memory for caching too... this is likely moot unless very large arrays are required
# https://joblib.readthedocs.io/en/latest/memory.html
# From the docs, this should be fast for large arrays, while other cache methods may not be.
# Q: is it possible to set cache to RAM?
# cachedir = r'E:\temp\memCache'  # SSD cache
# memory = Memory(cachedir, verbose=0)
# @memory.cache
# def Wigner3jMemory(QNlist):
#     """Wrapper for 3j caching with functools.lru_cache"""
#     return sf.Wigner3j(QNlist[0], QNlist[1], QNlist[2], QNlist[3], QNlist[4], QNlist[5])

# Parallelised 3j using Numba prange, see sf_function_tests_110220.py
@nb.njit(parallel = True)
def w3jprange(QNs):
    """
    Wrapper for 3j with @numba.njit(parallel=True), using prange parallel loop.

    In testing (Feb 2020) on an AMD Threadripper 1950X (16 core) this was (usually) fastest case, and beat vectorised version.

    Parameters
    ----------
    QNs : np.array
        Array of QNs to calculated Wigner 3j terms for, columns [l,lp,L,m,mp,M].

    Returns
    -------
    w3j_QNs : np.array
        Array of Wigner 3j results, one per row of input QNs.

    """

    w3j_QNs = np.zeros(QNs.shape[0])
    for n in nb.prange(0, QNs.shape[0]):
        w3j_QNs[n] = sf.Wigner3j(QNs[n,0], QNs[n,1], QNs[n,2], QNs[n,3], QNs[n,4], QNs[n,5])
    return w3j_QNs


# To try - other high-level parallelisation methods, e.g. xyzpy parallel, joblib

# Vectorized 3j, see sf_function_tests_110220.py
# Use Numba to compile function, vecotrised over input sets of QNs, doesn't return, but writes directly to passed 1D list of corresponding 3j terms
@nb.guvectorize(["void(int32[:,:], float64[:])"], '(n,m)->(n)', target = 'parallel')
def w3jguVecCPU(QNs,w3j_QNs):
    """
    Wrapper for 3j with vectorization via @numba.guvectorize(["void(int32[:,:], float64[:])"], '(n,m)->(n)', target = 'parallel').

    Parameters
    ----------
    QNs : np.array
        Array of QNs to calculated Wigner 3j terms for, columns [l,lp,L,m,mp,M].

    w3j_QNs : np.array
        Empty array to hold results (no return from @guvectorize).
        Create as w3j_QNs = np.zeros(QNs.shape[0])

    """
    for n in range(QNs.shape[0]):
        w3j_QNs[n] = sf.Wigner3j(QNs[n,0], QNs[n,1], QNs[n,2], QNs[n,3], QNs[n,4], QNs[n,5])


# # Try wrapping vec version for caching too...
# @lru_cache(maxsize = None)
# def Wigner3jguVecCPULRUCached(QNs,w3j_QNs):
#     w3jguVecCPU(QNs,w3j_QNs)
# # def Wigner3jguVecCPULRUCached(QNs):
# #     w3jguVecCPU(QNs)
#
# @memory.cache
# def Wigner3jguVecCPUMemory(QNs,w3j_QNs):
#     w3jguVecCPU(QNs,w3j_QNs)
# # def Wigner3jguVecCPUMemory(QNs):
# #     w3jguVecCPU(QNs,w3j_QNs)
