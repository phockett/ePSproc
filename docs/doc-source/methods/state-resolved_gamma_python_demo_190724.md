---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# State-resolved gamma calculations (python)
19/07/24

+++

For calculations, making use of state-to-state $\gamma$ and $C$ parameters, some basic routines are available in :py:module:`epsproc.geomCalc.gamma`.

- For basic use, working from pre-tabulated (legacy) gamma params, [see the "legacy" notes](../state-resolved_gamma_legacy_demo_170724.html).
- In this notebook, computation of gamma params with python routines is demonstrated.

The C and gamma parameters and computations are defined as per Refs. [1,2], in particular Ref. [2], Sect 3.1:

$$
\begin{eqnarray}
C(lm\lambda N_{t}M_{i}\mu_{\lambda}) & = & (2N_{t}+1)(-1)^{M_{+}+q}\left(\begin{array}{ccc}
N_{t} & 1 & l\\
M_{t} & p & m
\end{array}\right)\left(\begin{array}{ccc}
N_{+} & N_{i} & N_{t}\\
-M_{+} & M_{i} & M_{t}
\end{array}\right)\nonumber \\
 & \mathsf{x} & \left(\begin{array}{ccc}
N_{+} & N_{i} & N_{t}\\
-K_{+} & K_{i} & K_{t}
\end{array}\right)\left(\begin{array}{ccc}
N_{t} & 1 & l\\
-K_{t} & q & -\lambda
\end{array}\right)\nonumber \\
 & \mathsf{x} & \left(\begin{array}{ccc}
N_{+} & J_{+} & S_{+}\\
M_{+} & M_{J+} & M_{S+}
\end{array}\right)\left(\begin{array}{ccc}
N_{+} & J_{+} & S_{+}\\
K_{+} & P_{+} & \Sigma_{+}
\end{array}\right)\label{eq:geom-params-C}
\end{eqnarray}
$$

$$
\begin{eqnarray}
\gamma_{\alpha\alpha_{+}l\lambda ml'\lambda'm'} & = & (2N_{i}+1)(2N_{+}+1)(-i)^{l'-l}\sum_{M_{+}}\sum_{M_{i}M_{i}'}\sum_{N_{t}N_{t}'}\sum_{\mu_{\lambda}\mu_{\lambda}'}{}^{J_{i}K_{i}}\boldsymbol{\rho}_{M_{i}M_{i}'}\nonumber \\
 & \mathsf{x} & C(lm\lambda N_{t}M_{i}q)C(l'm'\lambda'N_{t}'M_{i}'q')\label{eq:gamma-state}
\end{eqnarray}
$$

This form allows for:

- Angular momentum transfer between initial ($_i$) and final ($_+$) states, with transfer terms denoted by subscript $t$.
- Modulation by initial $M$-state distribution, expressed as a density matrix $\boldsymbol{\rho}_{M_{i}M_{i}'}$.
- Cf. aligned-frame case, which assumes a _decoupled_ rotational wavepacket, and is derived from a sum over $J,M$ states (see Ref. [3] for further details).


Refs:

1. Hockett, Paul. 2009. “Photoionization Dynamics of Polyatomic Molecules.” PhD Thesis, University of Nottingham. http://eprints.nottingham.ac.uk/10857/.
2. ———. 2018. Quantum Metrology with Photoelectrons, Volume 1: Foundations. IOP Publishing. https://doi.org/10.1088/978-1-6817-4684-5.
3. Stolow, Albert, and Jonathan G. Underwood. 2008. “Time-Resolved Photoelectron Spectroscopy of Non-Adiabatic Dynamics in Polyatomic Molecules.” In Advances in Chemical Physics, edited by Stuart A. Rice, 139:497–584. Advances in Chemical Physics. Hoboken, NJ, USA: John Wiley & Sons, Inc. https://doi.org/10.1002/9780470259498.ch6.

+++

## Gamma calc demo

```{code-cell} ipython3
# Load functions
from epsproc.geomFunc.gamma import gammaCalc
```

```{code-cell} ipython3
# Compute for default case
gammaCalc.gammaCalc()
```

```{code-cell} ipython3
# Compute for specific channel = [Ni,Ki,Nc,Kc]

channel = [2,0,1,1]
gammaCalc.gammaCalc(channel)
```

```{code-cell} ipython3
# Compute for specific channel = [Ni,Ki,Nc,Kc]
# + include initial state density matrix
# NOTE: TODO, density matrix not yet implemented.

channel = [2,0,1,1]
gammaCalc.gammaCalc(channel, denMat=1)
```

## Versions

```{code-cell} ipython3
import scooby
scooby.Report(additional=['epsproc', 'holoviews', 'hvplot', 'xarray', 'matplotlib', 'bokeh'])
```

```{code-cell} ipython3
# Check current Git commit for local ePSproc version
import epsproc as ep
from pathlib import Path
!git -C {Path(ep.__file__).parent} branch
!git -C {Path(ep.__file__).parent} log --format="%H" -n 1
```

```{code-cell} ipython3
# Check current remote commits
!git ls-remote --heads https://github.com/phockett/ePSproc
```

```{code-cell} ipython3

```
