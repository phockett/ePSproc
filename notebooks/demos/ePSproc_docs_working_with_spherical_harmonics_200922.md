---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

+++ {"tags": []}

# Working with spherical harmonics
20/09/22

This notebook gives a brief introduction to some of the conventions and functionality for handling spherical harmonics and related functions, using low-level ePSproc routines. Note some higher-level routines are available in the [base class](../demos/ePSproc_class_demo_161020.html).

## Expansions in complex spherical harmonics 

In general, expansions in complex harmonics are used, with expansion parameters $\beta_{L,M}$. A given function is then written as:

\begin{equation}
I(\theta,\phi)=\sum_{l,m}\beta_{l,m}Y_{l,m}(\theta,\phi)\label{eq:sph-harmonics}
\end{equation}

The default routines in ePSproc [make use of scipy.special.sph_harm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html) as the default calculation routine. Note the $(\theta, \phi)$ definition, and normalisation, corresponding to the usual "physics" convention ($\theta$ defined from the z-axis, includes the Condon-Shortley phase, and orthonormalised):

\begin{equation}
Y_{l,m}(\theta,\phi) = (-1)^m\sqrt{\frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!}}e^{i m \phi} P^m_l(\cos(\theta))
\end{equation}


For more details, see [also Wikipaedia, which has matching defintions plus further discussion](https://en.wikipedia.org/wiki/Spherical_harmonics). Note that in the Scipy routines the Condon-Shortley phase term, $(-1)^m$, is actually included in the associated Legendre polynomial function $P^m_l$, but is written explicitly above.

For additional control and conversion, the [SHtools library can also be used](https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html) - this includes some useful utility functions including converting between different forms (e.g. real and complex forms), and basic plotters.

## Expansions in symmetrized harmonics

In some case, symmetrized or generalised harmonics are useful. These are defined similarly, but with expansion coeffs additionally set by (point group) symmetry. A [basic implementation can be found in the PEMtk package](https://pemtk.readthedocs.io/en/latest/sym/pemtk_symHarm_demo_160322_tidy.html), defined as:

\begin{equation}
X_{hl}^{\Gamma\mu*}(\theta,\phi)=\sum_{\lambda}b_{hl\lambda}^{\Gamma\mu}Y_{l,\lambda}(\theta,\phi)\label{eq:symm-harmonics}
\end{equation}

Where the $b_{hl\lambda}^{\Gamma\mu}$ are defined by symmetry. General function expansions can then be written as a set of spherical harmonics including symmetrization, or an equivalent expansion in symmetrized harmonics (here not all indices may be necessary):

\begin{equation}
I(\theta,\phi)=\sum_{h,l,\Gamma,\mu}\beta_{hl}^{\Gamma\mu}X_{hl}^{\Gamma\mu*}(\theta,\phi)=\sum_{h,l,\Gamma,\mu}\sum_{\lambda}b_{hl\lambda}^{\Gamma\mu}\beta_{l,\lambda}Y_{l,\lambda}(\theta,\phi)\label{eq:symm-harmonics-2}
\end{equation}

## Expansions in Legendre polynomials

For cylindrically-symmetric cases ($m=0$), general expansions in Legendre polynomials are also sufficient and often used:

\begin{equation}
I(\theta)=\sum_{l}\beta_{l}P^0_l(\cos(\theta))\label{eq:lg-expansion}
\end{equation}

Where the "Legendre" and "spherical harmonic" expansion parameters for the harmonics defined above can be related by:

\begin{equation}
\beta^{Sph}_{l,0} = \sqrt{(2L+1)/4\pi}\beta^{Lg}_{l}
\end{equation}

+++

## Numerical implementation

+++

### Imports

```{code-cell} ipython3
import xarray as xr
import pandas as pd
import numpy as np
import epsproc as ep
```

### Computing spherical harmonics on a grid

A basic wrapper for the backends is provided by `ep.sphCalc`. This computes spherical harmonics for all orders up to `Lmax`, for a given angular resolution or set of angles, and returns an Xarray.

```{code-cell} ipython3
# Compute harmonics for 50x50 (theta,phi) grid
Isph = ep.sphCalc(Lmax = 2, res = 50)
Isph
```

```{code-cell} ipython3
# Compute harmonics for 50x25 (theta,phi) grid
I = ep.sphCalc(Lmax = 2, res = [50,25])
I
```

```{code-cell} ipython3
# These can be plotted ... although may not be very useful
ep.sphSumPlotX(I)
```

```{code-cell} ipython3
# Change plot style to default MPL
ep.basicPlotters.setPlotters(resetMpl=True)
```

```{code-cell} ipython3
# For Legendre Polynomials, set fnType = 'lg'
Ilg = ep.sphCalc(Lmax = 2, res = 50, fnType = 'lg')
ep.sphSumPlotX(Ilg)
```

```{code-cell} ipython3
Ilg
```

### Setting coefficients

To manually set coefficients from lists or arrays, use the `setBLMs()` function.

```{code-cell} ipython3
# Set specific LM coeffs by list with setBLMs, items are [l,m,value]
from epsproc.sphCalc import setBLMs
BLM = setBLMs([[0,0,1],[1,1,1],[2,2,1]])   # Note different index
BLM
```

```{code-cell} ipython3
# Set specific LM coeffs by numpy array with setBLMs, items are [l,m,value]
from epsproc.sphCalc import setBLMs
BLM = setBLMs(np.array([[0,0,1],[1,1,1],[2,2,1]]))   # Note different index
BLM
```

```{code-cell} ipython3
# Set all allowed BLM up to Lmax with genLM
LM = ep.genLM(2)
# BLMall = setBLMs(np.ones([LM.shape[0],1]), LM)
BLMall = setBLMs(np.ones(LM.shape[0]), LM)
BLMall
```

```{code-cell} ipython3
# Set all allowed BLM up to Lmax with genLM, random value case
LM = ep.genLM(2)
BLMallrand = setBLMs(np.random.random(LM.shape[0]), LM)
BLMallrand
```

Note that the function adds a `t` dimension, which is used as a generic index in these cases. For 1D cases it can be dropped by squeezing the array:

```{code-cell} ipython3
BLMallrand = BLMallrand.squeeze(drop=True)
BLMallrand
```

... or, for cases $\beta_{L,M}(t)$, can be set directly (or renamed and used for other coords). Currently only a single dimension is supported here.

```{code-cell} ipython3
# Set all allowed BLM up to Lmax with genLM, random value case, and as a function of t
# If coordinates are not specified, integer labels are used.
LM = ep.genLM(2)
BLMallrandt = setBLMs(np.random.random([LM.shape[0],5]), LM)
BLMallrandt
```

### Calculating & plotting expansions

This can be done manually with Xarray multiplications for full control, or via the function `sphFromBLMPlot()`.

```{code-cell} ipython3
# Direct Xarray tensor multiplication - this will correctly handle dimensions if names match.
Irand = BLMallrand.rename({'BLM':'LM'}) * Isph
# ep.sphSumPlotX(Irand, facetDim = 't')   # Note this may need facetDim set explicitly here for non-1D cases
ep.sphSumPlotX(Irand)
Irand
```

```{code-cell} ipython3
# Functional form - this returns array, include plotFlag = True to show plot and return a fig handle
Itp, fig = ep.sphFromBLMPlot(BLMallrand, plotFlag = True)
Itp
```

```{code-cell} ipython3
# Set the backend to 'pl' for an interactive surface plot with Plotly
ep.sphFromBLMPlot(BLMallrandt, facetDim = 't', plotFlag = True, backend = 'pl');
```

## Conversions

+++

### Coefficient conversion

To convert $\beta_{L,M}$ normalisations and types, a few utility functions are implemented.

```{code-cell} ipython3
# Set some terms - note the BLM values are initially unlabelled by normType or with 'harmonics' attribututes.
LM = ep.genLM(Lmax =2, allM=False)
Blg = setBLMs(np.random.random(LM.shape[0]), LM, dtype = 'Lg', kind = 'real')  # Additional labels are set in Blg.attrs['harmonics']
#.sel(m=0, drop=False)   # Testing in XR 2022.3.0, this ALWAYS drops m for stacked dim case - should add some dim handling to conv_BL_BLM()!
Blg
```

```{code-cell} ipython3
# Convert from Sph to Lg normalised values and don't renormalise
Bsp = ep.util.conversion.conv_BL_BLM(Blg, to='sph', renorm=False)
Bsp
```

### SHtools conversion

To convert values to [SHtools objects](https://shtools.github.io/SHTOOLS/index.html), use `epsproc.sphFuncs.sphConv.SHcoeffsFromXR()`. Note this currently (Sept. 2022) only supports 1D input arrays.

```{code-cell} ipython3
BLMallrand
```

```{code-cell} ipython3
from epsproc.sphFuncs.sphConv import *

# Create SHtools object - note some additional settings may need to be specified
sh = SHcoeffsFromXR(BLMallrand)   #, kind = 'real')
sh
```

```{code-cell} ipython3
# Check parameters - this are set as a 3D array, where the SHtools convention is to  index as [+/-,l,m]
sh.coeffs
```

```{code-cell} ipython3
# SHtools has various conversion routines
sh.convert(csphase =1 , normalization='ortho')
```

```{code-cell} ipython3
# SHtools expand on grid & plot
lmaxGrid = 10
grid = sh.expand(lmax = lmaxGrid)
grid.plot(colorbar='right')   
```

## Versions

```{code-cell} ipython3
import scooby
scooby.Report(additional=['epsproc', 'holoviews', 'hvplot', 'xarray', 'matplotlib', 'bokeh'])
```

```{code-cell} ipython3
# Check current Git commit for local ePSproc version
from pathlib import Path
!git -C {Path(ep.__file__).parent} branch
!git -C {Path(ep.__file__).parent} log --format="%H" -n 1
```

```{code-cell} ipython3
# Check current remote commits
!git ls-remote --heads https://github.com/phockett/ePSproc
```
