"""
ePSproc electric field functions.

Generate, propagate and model E-fields.


06/05/20    Updated methods, added Xarray conversion & Holoviews plotting - basics only.

27/04/20    Updating & packaging - but still very much a work-in-progress.

20/3/20     Basics developed in tests/methodDev/geometric_method_dev_low-level_E-fields_200320.ipynb

TODO
- Integration with EPR and geometric methods code. Need Xarray conversion plus appropriate regridding for this, and proper handling of polarisation terms p.
- Tidy up and debug.
- Improve spectrogram methods & plotting.
- Improve pulse modification methods.
- Additional pulse properties, e.g. instantaneous frequency.


"""


# Tags for mkinit, https://github.com/Erotemic/mkinit
# On Stimpy, run with
# > mkinit D:\code\github\ePSproc\epsproc > mkinitOut.txt
# <AUTOGEN_INIT>
from epsproc.efield import efields

from epsproc.efield.efields import (Efield,)

__all__ = ['Efield', 'efields']
# </AUTOGEN_INIT>
