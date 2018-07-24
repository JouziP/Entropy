# -*- coding: utf-8 -*-
import numpy as np

def getEnergy(lattice, neighbs, weights, J, H):
    N1 = lattice.shape[0]
    N2 = lattice.shape[1]
    
    E = 0
    
    for i2 in range(N2):
        for i1 in range(N1):
            s = lattice[ i1, i2]
            s_idx  =  i2 * N1 + i1
#            print s_idx
            E += H * s * 2
            
            for k in range(1, neighbs.shape[0]):
                neighb_idx = int(neighbs[ k,  s_idx])
                weight = weights[k, s_idx]
                
                j2 = neighb_idx/N1
                j1 = neighb_idx - j2*N1
                
#                print s_idx, neighb_idx,  j1, j2, weight
                s_neighb = lattice[j1, j2]
                
                E += J * weight * s * s_neighb 
            
    E = E/2.  # double counting
    return E
    