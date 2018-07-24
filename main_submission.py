# -*- coding: utf-8 -*-

## -*- coding: utf-8 -*-
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from initialize import initialize
from getEnergy import getEnergy
from evolveInTemp import evolveInTemp
from getTotalPolarization import getTotalPolarization
from getEntropyFromHC import getEntropyFromHC
from buildLattice import buildLattice

#### inputs
## seed
SEED = 113
np.random.seed(SEED)
## lattice
N1 = int(sys.argv[1])
N2 = int(sys.argv[2])
neigb_cutOff = int(sys.argv[6])  # max rank of neighbor: 1 nn only, 2 nn and nnn only etc.
theta_factor = float(sys.argv[4])  #
a1 = 1
a2 = float(sys.argv[7])
w0 = 1
power = float(sys.argv[8])
# Hamilt
J = float(sys.argv[3])#-1 # exchange coupling
H = float(sys.argv[5]) # external field
# Temps
T_min = 0.001
T_max =  int( float(sys.argv[9])   *N1*N2) # 40
num_T = int(float(sys.argv[10])    * T_max) # 400
## MC
num_warmup = 1000
num_sampling = 200
num_mc = 100
########################

alpha = np.pi/ theta_factor
# build lattice 
w_clean, n_clean = buildLattice(N1, N2, alpha, neigb_cutOff)
lattice = np.zeros([N1, N2])





# initialize random
initialize(lattice)
# must use neighbs and weight now
E_int  = getEnergy(lattice, n_clean, w_clean, J, H)
P_int =  getTotalPolarization(lattice)
print '\n\n    ********   ',   E_int, P_int
print lattice
# temps

dT = (T_max - T_min)/(num_T + 0.0)
Temps = np.array([ T_min + t*dT for t in range(num_T) ])
Temps = list(Temps[::-1])
### for t in Temps ...
#Temp = Temps[0]

E_vs_T = np.zeros([len(Temps), 2])
E2_vs_T = np.zeros([len(Temps), 2])
P_vs_T = np.zeros([len(Temps), 2])
P2_vs_T = np.zeros([len(Temps), 2])
HC_T = np.zeros([len(Temps), 2])
Ssp_T = np.zeros([len(Temps), 2])

#
for Temp in Temps:
    ## warm-up
    E_int, P_int = evolveInTemp(lattice, n_clean, w_clean, Temp,J, H, E_int, P_int, num_warmup)
#    print '\n\n    ********   ',  Temp, E_int, P_int
#    print lattice
    ## sampling
    # energy
    energy_avg=[]
    energy_avg.append(E_int)
    # energy ^ 2
    energy2_avg=[]
    energy2_avg.append(E_int*E_int)
    # polarization
    polarization_avg=[]
    polarization_avg.append(P_int)    
    # polarization ^ 2
    polarization2_avg=[]
    polarization2_avg.append(P_int*P_int)      
#    
    ###########################################    
   
    ##    
    for s in range(num_sampling):
        E_int, P_int = evolveInTemp(lattice, n_clean, w_clean, Temp,J, H, E_int, P_int, num_mc)
        ## sampling
        # energy
        energy_avg.append(E_int)
        # energy ^ 2
        energy2_avg.append(E_int*E_int)
        # polarization
        polarization_avg.append(P_int)    
        # polarization ^ 2
        polarization2_avg.append(P_int*P_int)  
        
    E_avg = np.average(energy_avg)
    E_vs_T[Temps.index(Temp), 0] = Temp
    E_vs_T[Temps.index(Temp), 1] = E_avg /(N1*N2 + 0.0)
    
    E2_avg = np.average(energy2_avg)
    E2_vs_T[Temps.index(Temp), 0] = Temp
    E2_vs_T[Temps.index(Temp), 1] = E2_avg/((N1*N2 + 0.0)**2)
    
    P_avg = np.average(polarization_avg)
    P_vs_T[Temps.index(Temp), 0] = Temp
    P_vs_T[Temps.index(Temp), 1] = P_avg/(N1*N2 + 0.0)
    
    P2_avg = np.average(polarization2_avg)
    P2_vs_T[Temps.index(Temp), 0] = Temp
    P2_vs_T[Temps.index(Temp), 1] = P2_avg/((N1*N2 + 0.0)**2)
    
    HC = -(E_avg**2 - E2_avg)/(N1*N2 + 0.0)/(Temp**2)
    HC_T[Temps.index(Temp), 0] = Temp
    HC_T[Temps.index(Temp), 1] = HC
    
    Ssp= -(P_avg**2 - P2_avg)/(N1*N2 + 0.0)/(Temp)
    Ssp_T[Temps.index(Temp), 0] = Temp
    Ssp_T[Temps.index(Temp), 1] = Ssp
    
    
######################## entropy 
S_T  = np.zeros([len(Temps), 2])
entropy = getEntropyFromHC(HC_T)
S_T[:, 0] = Temps
S_T[:, 1] = entropy
print S_T[-10:-1]

 
#plt.plot(E_vs_T[:, 0], E_vs_T[:, 1], '-o')
##plt.plot(E2_vs_T[:, 0], E2_vs_T[:, 1], '-o')
#plt.plot(P_vs_T[:, 0], P_vs_T[:, 1], '-o')
#plt.plot(P2_vs_T[:, 0], P2_vs_T[:, 1], '-o')
#plt.plot(Ssp_T[:, 0],Ssp_T [:, 1], '-o')
#plt.plot(HC_T[:, 0],HC_T [:, 1], '-o')
#plt.plot(S_T[:, 0],S_T [:, 1], '-o')



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
                                           
np.savetxt('E_vs_T'+filename+'.txt', E_vs_T, fmt='%2.4f\t%2.4f')
np.savetxt('E2_vs_T'+filename+'.txt', E2_vs_T, fmt='%2.4f\t%2.4f')
np.savetxt('P_vs_T'+filename+'.txt', P_vs_T, fmt='%2.4f\t%2.4f')
np.savetxt('Ssp_T'+filename+'.txt', Ssp_T, fmt='%2.4f\t%2.4f')
np.savetxt('HC_T'+filename+'.txt', HC_T, fmt='%2.4f\t%2.4f')
np.savetxt('S_T'+filename+'.txt', S_T, fmt='%2.4f\t%2.4f')


















