import numpy as np

def getEntropyFromHC(HC_T):
    Temps = HC_T[:, 0]
    specific_heat = HC_T[:, 1]    
    
    Entropy_array=[]
    
    Entropy=np.log(2)
    Entropy_array.append(Entropy)
    
    for T in range(len(Temps)-1):
        dEntropy=0
        for t in range(T, 0 , -1):
            dT = Temps[t-1] - Temps[t] 
            dEntropy+= dT * specific_heat[t]/Temps[t]
        Entropy=   ( np.log(2) - dEntropy  )  #+ 1./Temps[0] * E0
        Entropy_array.append(Entropy)
    return Entropy_array
