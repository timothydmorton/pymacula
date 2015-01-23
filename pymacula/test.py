from __future__ import print_function, division

import numpy as np
import numpy.random as rand
from ._macula import maculamod as macula


def test_macula():
    """Trying to reproduce example in maculacall.f90
    """
    ndata = 1e3
    Nspot = 5
    mmax = 1
    derivatives = True
    temporal = True
    TdeltaV = True

    t = np.linspace(0.1,100,1000)
    Tstart = 0; Tend = 100.01

    theta_star = [np.pi/2, 0.233*(Tend-Tstart) + 0.1*(Tend-Tstart),
                  0.2, 0.2,
                  0.4, 0.4269, -0.0227, -0.0839,
                  0.4, 0.4269, -0.0227, -0.0839]

    theta_spot = np.zeros((8,Nspot))
    for k in range(Nspot):
        theta_spot[0,k] = rand.random()*np.pi # longitude
        theta_spot[1,k] = rand.random()*np.pi/2 # latitude
        theta_spot[2,k] = rand.random()*10*np.pi/180 # alpha_max (size of spot)
        theta_spot[3,k] = rand.random()*0.5 # fspot (flux contrast)
        theta_spot[4,k] = rand.random()*(Tend-Tstart) + Tstart # tmax
        theta_spot[5,k] = rand.random()*(Tend-Tstart)
        theta_spot[6,k] = rand.random()*(Tend-Tstart)
        theta_spot[7,k] = rand.random()*(Tend-Tstart)

    theta_inst = np.array([1.,1.])

    return macula.macula(t, derivatives, temporal, TdeltaV, theta_star, theta_spot, theta_inst, Tstart, Tend) # fmod, dfmod_star, dfmod_spot, dfmod_inst, dfmoddt, deltaratio
    
    
