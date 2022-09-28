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

# Working with real spherical harmonics
20/09/22

This notebook gives a brief introduction to some of the conventions and functionality for handling *real* spherical harmonics and related functions, using low-level ePSproc routines. For *complex* spherical harmonics, and more general method notes, see ['working with complex harmonics'](ePSproc_docs_working_with_spherical_harmonics_200922.html) - note the complex forms are generally used in ePSproc, and general correct handling for cases using real spherical harmonics is not currently implemented/guaranteed.

Note some higher-level routines are available in the [base class](../demos/ePSproc_class_demo_161020.html).

## Definitions

Per [the Wikipaedia definitions](https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form), the real spherical harmonics can be defined as:

\begin{aligned}
Y_{\ell m}&={\begin{cases}{\dfrac {i}{\sqrt {2}}}\left(Y_{\ell }^{m}-(-1)^{m}\,Y_{\ell }^{-m}\right)&{\text{if}}\ m\lt0
\\Y_{\ell }^{0}&{\text{if}}\ m=0
\\{\dfrac {1}{\sqrt {2}}}\left(Y_{\ell }^{-m}+(-1)^{m}\,Y_{\ell }^{m}\right)&{\text{if}}\ m\gt0.\end{cases}}
\\&={\begin{cases}{\dfrac {i}{\sqrt {2}}}\left(Y_{\ell }^{-|m|}-(-1)^{m}\,Y_{\ell }^{|m|}\right)&{\text{if}}\ m\lt0
\\Y_{\ell }^{0}&{\text{if}}\ m=0
\\{\dfrac {1}{\sqrt {2}}}\left(Y_{\ell }^{-|m|}+(-1)^{m}\,Y_{\ell }^{|m|}\right)&{\text{if}}\ m\gt0.\end{cases}}
\\&={\begin{cases}{\sqrt {2}}\,(-1)^{m}\,\Im [{Y_{\ell }^{|m|}}]&{\text{if}}\ m\lt0
\\Y_{\ell }^{0}&{\text{if}}\ m=0
\\{\sqrt {2}}\,(-1)^{m}\,\Re [{Y_{\ell }^{m}}]&{\text{if}}\ m\gt0.\end{cases}}
\end{aligned}

Where $Y_{\ell m}$ are real harmonics, and $Y_{\ell }^{m}$ are complex.

And the inverse relations:

\begin{aligned}
Y_{\ell }^{m}={\begin{cases}{\dfrac {1}{\sqrt {2}}}\left(Y_{\ell |m|}-iY_{\ell ,-|m|}\right)&{\text{if}}\ m\lt0
\\[4pt]Y_{\ell 0}&{\text{if}}\ m=0
\\[4pt]{\dfrac {(-1)^{m}}{\sqrt {2}}}\left(Y_{\ell |m|}+iY_{\ell ,-|m|}\right)&{\text{if}}\ m\gt0.\end{cases}}
\end{aligned}

The real harmonics can also be written explicitly as:

\begin{aligned}
Y_{\ell m}={\begin{cases}\left(-1\right)^{m}{\sqrt {2}}{\sqrt {{\dfrac {2\ell +1}{4\pi }}{\dfrac {(\ell -|m|)!}{(\ell +|m|)!}}}}\;P_{\ell }^{|m|}(\cos \theta )\ \sin(|m|\varphi )&{\text{if }}m\lt0
\\[4pt]{\sqrt {\dfrac {2\ell +1}{4\pi }}}\ P_{\ell }^{m}(\cos \theta )&{\text{if }}m=0
\\[4pt]\left(-1\right)^{m}{\sqrt {2}}{\sqrt {{\dfrac {2\ell +1}{4\pi }}{\dfrac {(\ell -m)!}{(\ell +m)!}}}}\;P_{\ell }^{m}(\cos \theta )\ \cos(m\varphi )&{\text{if }}m\gt0\,.\end{cases}}
\end{aligned}

For futher information, see [wikipeadia]((https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form)) or the [SHtools introduction to real spherical harmonics](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html).

## Numerical implementation

ePSproc currently (Sept. 2022) implements: 

- Generation of real harmonics.
- Conversion of real harmonic coefficients to complex harmonic coefficients.

General handling for real spherical harmonics is not guaranteed however, although - as per above - they can be defined simply by suitable combinations of complex harmonics. Note that the default routines in ePSproc [make use of scipy.special.sph_harm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html) as the default calculation routine, this is implemented with a default $(\theta, \phi)$ definition, and normalisation, corresponding to the usual "physics" convention ($\theta$ defined from the z-axis, includes the Condon-Shortley phase, and orthonormalised). See ['working with complex harmonics'](ePSproc_docs_working_with_spherical_harmonics_200922.html) for further discussion; the real harmonics are then defined from the complex forms following the final definition given above, where the Condon-Shortley phase is already included in the $Y_{\ell }^{|m|}$ terms:

\begin{aligned}
Y_{\ell m}&={\begin{cases}{\sqrt {2}}\,\Im [{Y_{\ell }^{|m|}}]&{\text{if}}\ m\lt0
\\Y_{\ell }^{0}&{\text{if}}\ m=0
\\{\sqrt {2}}\,\Re [{Y_{\ell }^{m}}]&{\text{if}}\ m\gt0.\end{cases}}
\end{aligned}

+++

## Imports

```{code-cell} ipython3
import xarray as xr
import pandas as pd
import numpy as np
import epsproc as ep

# Set compact XR repr
xr.set_options(display_expand_data = False)
```

## Computing spherical harmonics on a grid

A basic wrapper for the backends is provided by `ep.sphCalc`. This computes spherical harmonics for all orders up to `Lmax`, for a given angular resolution or set of angles, and returns an Xarray. Set `fnType = 'real'` for real harmonics.

```{code-cell} ipython3
# Compute harmonics for 50x50 (theta,phi) grid
Ire = ep.sphCalc(Lmax = 2, res = 50, fnType = 'real')
Ire
```

```{code-cell} ipython3
# Plot some components
# Ire = ep.sphCalc(Lmax = 2, res = 50, fnType = 'real')
ep.sphSumPlotX(Ire.sel(l=2), facetDim = 'm', backend = 'pl');
```

```{code-cell} ipython3
# Compare with complex case (default)
IC = ep.sphCalc(Lmax = 2, res = 50)
ep.sphSumPlotX(IC.sel(l=2), facetDim = 'm', backend = 'pl', titleString='Abs');
ep.sphSumPlotX(IC.sel(l=2), facetDim = 'm', backend = 'pl', pType = 'r', titleString='Real component');
ep.sphSumPlotX(IC.sel(l=2), facetDim = 'm', backend = 'pl', pType = 'i', titleString='Imaginary component');
```

Note that the complex harmonics have rotations between the real and imaginary components, and produce a cylindrically-symmetric absolute valued map (here computed as np.abs(I) $=\sqrt(a^2+b^2)$ for the complex case $Z=a+ib$), whilst the real harmonics have rotations between + and - m pairs. 

In general this can be viewed as a phase rotation in the complex plane; this can be visualised more directly by looking at a phase map, which - [by definition for the individual components](ePSproc_docs_working_with_spherical_harmonics_200922.html#Expansions-in-complex-spherical-harmonics) - is given by the term $e^{im\phi}$. (For summations over sets of harmonics this becomes more interesting!)

```{code-cell} ipython3
# Visualise the phase
ep.plotTypeSelector(IC.unstack('LM'), pType = 'phase').plot(y='Theta', x='Phi', row='l', col='m', robust=True)
```

## Setting coefficients

To manually set coefficients from lists or arrays, use the `setBLMs()` function, including `kind='real'` to specify an expansion in real harmonics. These can be used and expanded directly, or via the `sphFromBLMPlot()` function.

```{code-cell} ipython3
# Set all allowed BLM up to Lmax with genLM, random value case
from epsproc.sphCalc import setBLMs
LM = ep.genLM(2)
BLMallrand = setBLMs(np.random.random(LM.shape[0]), LM, kind = 'real')
BLMallrand
```

```{code-cell} ipython3
# Direct Xarray tensor multiplication - this will correctly handle dimensions if names match.
Irand = BLMallrand.rename({'BLM':'LM'}).squeeze(drop=True) * Ire
# ep.sphSumPlotX(Irand, facetDim = 't')   # Note this may need facetDim set explicitly here for non-1D cases
ep.sphSumPlotX(Irand, facetDim=None, backend='pl');
```

```{code-cell} ipython3
# Functional form - this returns array, include plotFlag = True to show plot and return a fig handle
Itp, fig = ep.sphFromBLMPlot(BLMallrand.squeeze(), plotFlag = True, backend='pl')
Itp
```

## Converting coefficients real to complex form

Basic conversion is implemented in `epsproc.sphFuncs.sphConv.sphRealConvert()`.

This currently (Sept. 2022) has two methods implemented for conversion:

- 'sh': Set terms as per SHtools definitions, https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
        Note this form has purely real or purely imaginary terms.
- 'std': Set +/-m terms from |m| real harmonics as per wikipedia defn, https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
        Note: this form currently assumes only +M or -M for any given term, and also sets conjugates upon conversion if incConj=True.
        This may be phase-rotated relative to SHtools convention.
        NOTE: THIS METHOD CURRENTLY FAILS - phase issues somewhere, produces incorrect expansions in complex harmonics.
        
A basic demo is given below. In general, plotting the original and converted forms is a good test of consistency here, and SHtools can also be used directly as an independent method.

```{code-cell} ipython3
# Basic distribution
BLM = setBLMs([[0,0,1],[1,1,1],[2,2,1]], kind='real').squeeze(drop=True)
BLM
```

```{code-cell} ipython3
backend = 'pl'
ItpRe, fig = ep.sphFromBLMPlot(BLM, plotFlag = True, backend = backend)   #, fnType = 'real')
```

```{code-cell} ipython3
from epsproc.sphFuncs.sphConv import *

BLMconvSH = sphRealConvert(BLM, method='sh')   #, addCSphase = False)
ItpCSH, fig = ep.sphFromBLMPlot(BLMconvSH, plotFlag = True, backend = backend, pType='a')
BLMconvSH
```

### SHtools conversion & plots

```{code-cell} ipython3
# Check vs. shtools routine
# Create SHtools object - note some additional settings may need to be specified
sh = SHcoeffsFromXR(BLM)   #, kind = 'real')
shC = sh.convert(kind = 'complex')   #, csphase =1 , normalization='4pi')  # Optional conversion options, will be taken from BLM.attrs['harmonics'] if not set

# shC.coeffs   # Display array
tabulateLM(XRcoeffsFromSH(shC))  # Display table
```

```{code-cell} ipython3
# Compare params directly - looks good
tabulateLM(BLMconvSH)
```

```{code-cell} ipython3
# Plot maps with SHtools

# SHtools expand on grid & plot
lmaxGrid = 10
grid = sh.expand(lmax = lmaxGrid)
grid.plot(colorbar='right')
```

```{code-cell} ipython3
grid = shC.expand(lmax = lmaxGrid)
grid.plot(colorbar='right')
```

```{code-cell} ipython3
# Test vs. ePSproc conversion
shC2 = SHcoeffsFromXR(BLMconvSH)
grid = shC2.expand(lmax = lmaxGrid)
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
