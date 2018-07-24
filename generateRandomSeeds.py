# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import numpy as np
#
# numRealizations
#
num_realizations = 1000

## seed
seed_master = 1021
np.random.seed(seed_master)
SEEDs=[]

for i in range(num_realizations):
    SEEDs.append(np.random.randint(10, 10000) )
    
np.savetxt('randomSeeds_masterSEED_%d.txt'%seed_master, 
           np.array(SEEDs))