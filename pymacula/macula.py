from __future__ import print_function, division

import numpy as np
import numpy.random as rand
from ._macula import maculamod


class Star(object):
    def __init__(self, incl=np.pi/2, Peq=30.,
                 kappa2=0.3, kappa4=0.3, 
                 c1=0.3999, c2=0.4269,
                 c3=-0.0227, c4=-0.839,
                 d1=0.3999, d2=0.4269,
                 d3=-0.0227, d4=-0.839):
        self.incl = incl
        self.Peq = Peq
        self.kappa2 = kappa2
        self.kappa4 = kappa4
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.d1 = d1
        self.d2 = d2
        self.d3 = d3
        self.d4 = d4
        

    @property
    def pars(self):
        return np.array([self.incl, self.Peq, self.kappa2, 
                         self.kappa4, self.c1, self.c2,
                         self.c3, self.c4, self.d1,
                         self.d2, self.d3, self.d4])
                 
class Spot(object):
    def __init__(self, lon=None, lat=None,
                 alpha_max=10., contrast=0.3,
                 tmax=None, lifetime=None,
                 ingress=None, egress=None):

        #assign random longitudes/latitudes if not provided
        if lon is None:
            lon = rand.random()*2*np.pi - np.pi
        self.lon = lon

        if lat is None:
            lat = np.arccos(rand.random())

        self.lat = lat

        self.alpha_max = np.deg2rad(alpha_max)
        self.contrast = contrast #fspot

        #all times are between 0 and 1; should be normalized
        # to actual data time span
        self.tmax = rand.random()
        self.lifetime = 1
        self.ingress = rand.random()
        self.egress = rand.random()

    @property
    def pars(self):
        return np.array([self.lon, self.lat,
                         self.alpha_max, self.contrast,
                         self.tmax, self.lifetime,
                         self.ingress, self.egress])

class MaculaModel(object):
    def __init__(self, t=None, fobs=None, 
                 tmin=0, tmax=500,
                 t_start=None, t_end=None,
                 star=None, nspots=3, spots=None,
                 inst_offsets=None, blend_factors=None):

        if np.size(t_start)!=np.size(t_end):
            raise ValueError('t_start and t_end must be same length')

        #number of data sets
        self.n_datasets = np.size(t_start)

        if inst_offsets is None:
            self.inst_offsets = np.ones(self.n_datasets)
        if blend_factors is None:
            self.blend_factors = np.ones(self.n_datasets)
            
        if t is None:
            self.t_start = np.array([tmin-0.01])
            self.t_end = np.array([tmax+0.01])
        else:
            if t_start is None:
                self.t_start = np.array([t[0] - 0.01])
            else:
                self.t_start = t_start

            if t_end is None:
                self.t_end = np.array([t[-1] + 0.01])
            else:
                self.t_end = t_end
        
        self.t_span = self.t_end[-1] - self.t_start[0]

        self.t = t
        self.fobs = fobs
        
        if star is not None:
            self.star = star
        else:
            self.star = Star()

        if spots is not None:
            self.spots = spots
        else:
            self.spots = [Spot() for i in xrange(nspots)]

        self.nspots = len(self.spots)

    def __call__(self, t, derivatives=False,
                 temporal=False, tdeltav=False,
                 full_output=False):
        """Calls macula

        Only works for t in defined range.
        """
        return macula(t, self.theta_star, self.theta_spot, self.theta_inst,
                      derivatives=derivatives, temporal=temporal,
                      tdeltav=tdeltav, tstart=self.t_start,
                      tend=self.t_end, full_output=full_output)
                      

    @property
    def theta_star(self):
        return self.star.pars

    @property
    def theta_spot(self):
        theta = np.array([self.spots[i].pars for i in xrange(self.nspots)]).T

        #rescale spot times from (0,1) to full data span
        theta[4:, :] *= self.t_span

        return theta

    @property
    def theta_inst(self):
        return np.array([self.inst_offsets, self.blend_factors])


def macula(t, theta_star, theta_spot, theta_inst,
           derivatives=False, temporal=False, tdeltav=False,
           full_output=False, tstart=None, tend=None):
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

    if tstart is None:
        tstart = t[0] - 0.01
    if tend is None:
        tend = t[-1] + 0.01

    res = maculamod.macula(t, derivatives, temporal, tdeltav, theta_star, theta_spot, theta_inst, tstart, tend) # fmod, dfmod_star, dfmod_spot, dfmod_inst, dfmoddt, deltaratio
    if full_output:
        return res
    else:
        return res[0] 
