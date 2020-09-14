# Verified CG term phase conventions
# These give correct B2 for N2 test case.

# l,m terms
phaseCons['lfblmCGCons']['negmp'] = False # mp > -mp  sign flip
                # Including this restricts things to mp=0 only? Phase con choice? Doesn't seem correct!
                # Full calculation including this term sends PU continuum to zero/NaN.
phaseCons['lfblmCGCons']['negM'] = True  # M > -M sign flip.

# Photon terms
phaseCons['lfblmCGCons']['negmup'] = False # mup > -mup sign flip.
phaseCons['lfblmCGCons']['negMP'] = False  # M > -M sign flip, photon term.
                                        # Originally thought this was a typo, but do need to set False as per manuscript!
