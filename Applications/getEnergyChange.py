# -*- coding: utf-8 -*-
import numpy as np


def getEnergyChange(lattice, n_clean, w_clean, i1,i2, J, H):
    N1 = lattice.shape[0]
    #N2 = lattice.shape[1]
    
    s = lattice[i1, i2]
    s_idx  =  i2 * N1 + i1  
    E_local = H * s
    
    for k in range(1, n_clean.shape[0]):
        neighb_idx = int(n_clean[ k,  s_idx])
        w  = w_clean[k, s_idx]
        
        j2 = neighb_idx/N1
        j1 = neighb_idx - j2*N1
        s_n = lattice[j1, j2 ]
    
        E_local +=   J * w * s * s_n 
        
    return -1 * 2 * E_local



#######################
  
##### test if all the neighbs are correct    
if __name__=='__main__':
    from getTotalPolarization import getTotalPolarization
    from initialize import initialize
    from getEnergy import getEnergy
    from buildLattice import buildLattice
    SEED = 113
    np.random.seed(SEED)
    
    N1 = 3
    N2 = 3
    neigb_cutOff = 1  # max rank of neighbor: 1 nn only, 2 nn and nnn only etc.
    theta_factor = 2  #
    alpha = np.pi/ theta_factor
    # build lattice 
    w_clean, n_clean = buildLattice(N1, N2, alpha, neigb_cutOff)
    lattice = np.zeros([N1, N2])
    
    
    J = 1 # exchange coupling
    H = 0  # external field
    
    
    # initialize random
    initialize(lattice)
    # must use neighbs and weight now
    E_int  = getEnergy(lattice, n_clean, w_clean, J, H)
    P_int =  getTotalPolarization(lattice)
    print '\n\n    ********   ',   E_int, P_int
    print lattice
    i1 = 0 
    i2 = 0
    print getEnergyChange(lattice, n_clean, w_clean, i1,i2, J, H)
    lattice[i1, i2] = - lattice[i1, i2]
    E_int  = getEnergy(lattice, n_clean, w_clean, J, H)
    P_int =  getTotalPolarization(lattice)
    print '\n\n    ********   ',   E_int, P_int
    print lattice