# -*- coding: utf-8 -*-

import numpy as np

def initialize(lattice):
    N1 = lattice.shape[0]
    N2 = lattice.shape[1]
    
    for i1 in range(N1):
        for i2 in range(N2):
            lattice[i1,i2] = (-1) ** np.random.randint(2)