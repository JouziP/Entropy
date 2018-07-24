# -*- coding: utf-8 -*-

import numpy as np
from Applications.getEnergyChange import getEnergyChange

def evolveInTemp(lattice, n_clean, w_clean, Temp, J, H, 
                  num_moves,
                 E_int, P_int, num_mc = 300):   
    
    N1 = lattice.shape[0]
    N2 = lattice.shape[1]
    
    for mc in range(num_mc):
        i1 = np.random.randint(N1)
        i2 = np.random.randint(N2)
        
        dE = getEnergyChange(lattice, n_clean, w_clean, i1, i2, J, H)
        E_fin  = E_int + dE
        
        p_trans = np.exp( -(E_fin - E_int)/(Temp+0.0) )
        r  = np.random.uniform(0,1)
#        print r, p_trans
        if p_trans > r:
            lattice[i1, i2] = - lattice[i1, i2]
            E_int = E_fin
            P_int += 2 * (lattice[i1, i2])
            num_moves+=1
        
    return E_int, P_int, num_moves
                