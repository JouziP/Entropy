# -*- coding: utf-8 -*-


def getTotalPolarization(lattice):
    N1 = lattice.shape[0]
    N2 = lattice.shape[1]
    
    polarization = 0

    for i1 in range(N1):
        for i2 in range(N2):
            polarization += lattice[i1,i2]
            
    return polarization