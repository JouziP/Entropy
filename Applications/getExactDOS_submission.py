# -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np

import sys
from Applications.getEnergy import getEnergy
from buildLattice import buildLattice



#################################
def updateLattice(lattice, config):
    N1 = lattice.shape[0]
    N2 = lattice.shape[1]
   
    for i2 in range(N2):
        for i1 in range(N1):
            s_idx  =  i2 * N1 + i1
            lattice[i1, i2] = (-1)**config[s_idx]
    
    
#################################    
def getBinaryArray(vecLen, num):
    maxNum = 2**vecLen-1
    if num<=maxNum:
        binaryVec = []
        idx=vecLen-1
        while (num >= 2):
            m= num/2
            r = num - 2*m
            binaryVec.append(int(r));
            num = m;
            idx-=1 
        binaryVec.append(int(num));
        while len(binaryVec)<vecLen:
            binaryVec.append(int(0))
        binaryVec.reverse()
        return binaryVec
    else:
        return    
    


#
# lattice
#
N1 = int(sys.argv[1])
N2 = N1
neigb_cutOff = int(sys.argv[2]) # max rank of neighbor: 1 nn only, 2 nn and nnn only etc.
theta_factor = float(sys.argv[3]) #
alpha = np.pi/theta_factor
a1 = 1
a2 = 1
w0 = 1
power = 3
#
# Hamilt
#
J = 1 # exchange/dipolar coupling
H = sys.argv[4]  # external field
#
# Temps
# for 10 -> 0.2, 10; 15 -> 
T_min = 0.001
T_max = int(0.2*N1*N2) # 40
num_T = int(10 * T_max) # 400
#
# MC
#
num_warmup = 300
num_sampling = 200
num_mc = 100
########################
# build lattice 
w_clean, n_clean = buildLattice(N1, N2, alpha, neigb_cutOff,
                                                        a1, 
                                                        a2,                  
                                                        w0, 
                                                        power)
lattice = np.zeros([N1, N2])

i =2**((N1*N2))-10
config=getBinaryArray(N1*N2, i) 

E_s = []
for i in range(2**(N1*N2)):
    config=getBinaryArray(N1*N2, i) 
    updateLattice(lattice, config)
    E = getEnergy(lattice, n_clean, w_clean, J, H)
    E_s.append(E)
#    print config, E
#    print np.min(E_s), np.max(E_s)
E_s = np.array(E_s).round(3)
unique_Es =np.array(list(set(list(E_s))))
unique_Es.sort()

rslt = np.zeros([len(unique_Es), 2])
    

pde = []
for E in unique_Es:
    num_E=len(np.where(E_s==E)[0])
    pde.append(num_E)
    
#plt.bar(unique_Es, pde)
rslt[:, 0]= unique_Es[:]
rslt[:, 1]= np.array(pde)[:]
    ######################################  save on file
#
filename = '_N1_%d_N2_%d_neigbCutOff_%d_thetaFactor_%1.1f_'%(N1, 
                                                           N2, 
                                                           neigb_cutOff, 
                                                           theta_factor)+\
        'a1_%d_a2_%d_w0_%1.1f_power_%1.1f'%(a1, a2, w0, power)+\
        'J_%1.1f_H_%1.1f_Tmin_%d_Tmax_%d_numT_%d_'%(J, 
                                                    H, 
                                                    T_min, 
                                                    T_max, 
                                                    num_T)+\
        'warmUp_%d_sampling_%d_mc_%d'%(num_warmup,
                                       num_sampling,
                                       num_mc)

np.savetxt('dos'+filename+'.txt', rslt, fmt='%2.4f\t%2.4f')


    
    
    
    
    
    
    
    
    
    
    
    
    