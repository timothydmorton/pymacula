from __future__ import print_function, division

import numpy as np
import numpy.random as rand
from ._macula import maculamod as macula


def MaculaModel(object):
    def __init__(self, t, fobs):
        self.t = t
        self.fobs = fobs

    

def macula(t, theta_star, theta_spot, theta_inst,
           derivatives=False, temporal=False, tdeltav=False,
           full_output=False):
    """Wrapper for macula FORTRAN routine.

    Parameters
    ----------
    theta_star : array_like
        Array of 12 parameters describing star:
        ! ------------------------------------------------------------------------------
        ! Theta_star(j) = Parameter vector for the star's intrinsic parameters
        ! ------------------------------------------------------------------------------
        ! Istar 	= Theta_star(1)		! Inclination of the star [rads]
        ! Peq 		= Theta_star(2)		! Rot'n period of the star's equator [d]
        ! kappa2 	= Theta_star(3)		! Quadratic differential rotation coeff
        ! kappa4 	= Theta_star(4)		! Quartic differential rotation coeff
        ! c1 		= Theta_star(5)		! 1st of four-coeff stellar LD terms
        ! c2 		= Theta_star(6)		! 2nd of four-coeff stellar LD terms
        ! c3 		= Theta_star(7)		! 3rd of four-coeff stellar LD terms
        ! c4 		= Theta_star(8)		! 4th of four-coeff stellar LD terms
        ! d1 		= Theta_star(9)		! 1st of four-coeff spot LD terms
        ! d2	 	= Theta_star(10)	! 2nd of four-coeff spot LD terms
        ! d3	 	= Theta_star(11)	! 3rd of four-coeff spot LD terms
        ! d4	 	= Theta_star(12)	! 4th of four-coeff spot LD terms
        
    theta_spot : array_like
        Array of spot parameters, shape (8, Nspot)
        ! ------------------------------------------------------------------------------
        ! Theta_spot(j,k) = Parameters of the k^th spot
        ! ------------------------------------------------------------------------------
        ! Lambda0(k) 	= Theta_spot(1,k)	! Longitude of spot at time tref(k)
        ! Phi0(k) 	= Theta_spot(2,k)	! Latitude of spot at time tref(k)
        ! alphamax(k)	= Theta_spot(3,k)	! Angular spot size at time tmax(k)
        ! fspot(k)	= Theta_spot(4,k)	! Spot-to-star flux contrast of spot k
        ! tmax(k)	= Theta_spot(5,k)	! Time at which spot k is largest
        ! life(k)	= Theta_spot(6,k)	! Lifetime of spot k (FWFM) [days]
        ! ingress(k)	= Theta_spot(7,k)	! Ingress duration of spot k [days]
        ! egress(k)	= Theta_spot(8,k)	! Egress duration of spot k  [days]

    theta_inst : array_like
        Nuisance/instrumental parameters for each of 'm' data sets.
        ! ------------------------------------------------------------------------------
        ! Theta_inst(j,m) = Instrumental/nuisance parameters
        ! ------------------------------------------------------------------------------
        ! U(m) 		= Theta_inst(1,m)	! Baseline flux level for m^th data set
        ! B(m) 		= Theta_inst(2,m)	! Blend factor for m^th data set

    derivatives : bool (optional)
        Whether to calculate derivatives.
    temporal : bool (optional)
        Whether to calculate temporal derivatives
    tdeltav : bool (optional)
        Whether to calculate transit depth variations
    full_output : bool (optional)
        If True, then return all output; otherwise just return model flux.
        
    """


    res = macula.macula(t, derivatives, temporal, TdeltaV, theta_star, theta_spot, theta_inst, Tstart, Tend) # fmod, dfmod_star, dfmod_spot, dfmod_inst, dfmoddt, deltaratio
    if full_output:
        return res
    else:
        return res[0] 
