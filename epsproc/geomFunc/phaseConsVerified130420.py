# Phase convention settings, case 'E' as verified for N2 (z,x,y) pol April 2020.
# Results as per https://epsproc.readthedocs.io/en/dev/methods/geometric_method_dev_pt2_170320_v060420.html

def setPhaseConventions(phaseConvention = 'S', typeList = False):
    """
    Set phase convention/choices for geometric functions.

    20/03/20 - first attempt. Aim to centralise all phase choices here to keep things clean and easy to debug/change.

    Set as dictionary for each term, to be appended to results Xarray.


    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, may have additional changes in numerics, e.g. conjugate Wigner D.
        If a dict of phaseConventions is passed they will simply be returned - this is for transparency/consistency over multiple fns which call setPhaseConventions()... although may be an issue in some cases.

    typeList : optional, bool, default = False
        If true, return list of supported options instead of list of phase choices.


    Note
    -----
    If a dict of phaseConventions is passed they will simply be returned - this is for transparency/consistency over multiple fns which call setPhaseConventions()... although may be an issue in some cases.

    """

    # Return just typeList if option set - this also defines master list.
    if typeList:
        # Supported types
        typeList = ['S', 'R', 'E']
        return typeList

    # If phaseConventions are preset, just return them.
    if type(phaseConvention) is dict:
        return phaseConvention


    # Set master dict to hold choices.
    phaseCons = {'phaseConvention':phaseConvention}

    #**** For generating QNs with genllpMatE()
    # Set here to avoid issues with dropped/missing terms later!
    # Some conventions will be tied to other choices below.
    genMatEcons = {}
    if phaseConvention == 'S':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp, otherwise M = -(m+mp)
        # genMatEcons['negM'] = False     # Set -M

    elif phaseConvention == 'R':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp, otherwise M = -(m+mp)
        # genMatEcons['negM'] = False     # Set -M

    elif phaseConvention == 'E':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp (normal convention), otherwise M = -(m+mp)
                                        # Note this is correlated with M switch terms later, incorrect settings may remove m or m' non-zero terms!

        # genMatEcons['negM'] = False     # Set -M

    phaseCons['genMatEcons'] = genMatEcons

    #*** For EPR tensor
    EPRcons = {}
    if phaseConvention == 'S':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'R':
        EPRcons['Rphase'] = True       # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'E':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    phaseCons['EPR'] = EPRcons

    #*** For Lambda term (as set by MFproj())
    lambdaCons = {}
    if phaseConvention == 'S':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'R':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'E':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = True  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    phaseCons['lambdaCons'] = lambdaCons


    #*** For Beta term (as set by betaTerm())
    betaCons = {}
    if phaseConvention == 'S':
        betaCons['negM'] = False       # Use -M term in 3j?
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    elif phaseConvention == 'R':
        betaCons['negM'] = False       # Use -M term in 3j?
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    elif phaseConvention == 'E':
        betaCons['negM'] = genMatEcons['negm']       # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    phaseCons['betaCons'] = betaCons


    #*** For MFPAD product case, as calculated in mfblmXprod()
    mfblmCons = {}
    if phaseConvention == 'S':
        mfblmCons['negRcoordSwap'] = True       # Swap -R and +R in EPR Xarray coords?

        mfblmCons['negMcoordSwap'] = True       # Swap +/-M coords.
        mfblmCons['Mphase'] = True              # Apply (-1)^M phase term.

        mfblmCons['negmCoordSwap'] = True       # Swap +/-m coords.
        mfblmCons['mPhase'] = True              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term.

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    if phaseConvention == 'R':
        mfblmCons['negRcoordSwap'] = True       # Swap -R and +R in EPR Xarray coords?

        mfblmCons['negMcoordSwap'] = True       # Swap +/-M coords.
        mfblmCons['Mphase'] = False              # Apply (-1)^M phase term.

        mfblmCons['negmCoordSwap'] = True       # Swap +/-m coords.
        mfblmCons['mPhase'] = True              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = False            # Apply (-1)^(mup - p) phase term. Already incoporated into MFproj()?

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    if phaseConvention == 'E':
        mfblmCons['negRcoordSwap'] = False       # Swap -R and +R in EPR Xarray coords? Already in EPRCons()?

        mfblmCons['negMcoordSwap'] = False       # Swap +/-M coords.
        mfblmCons['Mphase'] = betaCons['negM']   # Apply (-1)^M phase term. Correlated with +/-M term switch.

        mfblmCons['negmCoordSwap'] = False       # Swap +/-m coords.
        mfblmCons['mPhase'] = False              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term. Already incoporated into MFproj()?

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    phaseCons['mfblmCons'] = mfblmCons

    return phaseCons
