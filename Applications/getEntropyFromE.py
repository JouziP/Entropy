## external
import numpy as np

# internal
from getEntropyFromHC import getEntropyFromHC

def getEntropyFromE(E_vs_T):
    curve_smooth=getSmoothCurve(E_vs_T)   ### get a smooth curve
    HC_from_dervE_T = getDerivative(curve_smooth) # take derivative to obtain heat capacity
    
    Entropy_array = getEntropyFromHC(HC_from_dervE_T) # get entropy
    S_T=np.zeros([len(Entropy_array), 2])

    S_T[:, 0] = HC_from_dervE_T[:, 0]
    S_T[:, 1] = Entropy_array
    return S_T


def getDerivative(E_T):
    Temps = E_T[:, 0]
    E_derv_T = np.zeros([len(Temps)-1, 2])
    for t in range(len(Temps)-1):
        Temp = Temps[t]
        Temp_next = Temps[t+1]        
        dT = -1* (Temp_next - Temp)
        
        E = E_T[t, 1]
        E_next = E_T[t+1, 1]
        dE = -1* (E_next - E)
        
        E_derv = dE/(dT+0.0)
        
        E_derv_T[t, 0]= Temp
        E_derv_T[t, 1]= E_derv        
        
    return E_derv_T
        
def getSmoothCurve(curve_rough):
    Temps = curve_rough[:,  0]
    vals = curve_rough[: ,   1]
    
    intervals= 1
    curve_smoothed=np.zeros([len(Temps)-intervals, 2])
    curve_smoothed[0,1]=np.average(curve_rough[:intervals,1 ])
    curve_smoothed[0,0]=Temps[0]
    for t in range(intervals, len(Temps)-1):
        curve_smoothed[t-intervals+1,1] =\
                curve_smoothed[t-intervals,1] \
                +curve_rough[t, 1]/(intervals+0.0) \
                -curve_rough[t-intervals, 1]/(intervals+0.0)
        curve_smoothed[t-intervals+1,0]=Temps[t-intervals+1]
    return curve_smoothed
    
    
    
if __name__=='__main__':
    
    curve_smooth=getSmoothCurve(E_vs_T)
#    plt.plot(E_vs_T[:, 0], E_vs_T[:, 1], '-o')
#    plt.plot(curve_smooth[:, 0], curve_smooth[:, 1], '--')
    q= getDerivative(curve_smooth)
#    plt.plot(q[:, 0], q[:, 1], '--')
    S = getEntropyFromHC(q)
    
    S_T=np.zeros([len(S), 2])
    S_T[:, 1] = S
    S_T[:, 0] = q[:, 0]
    plt.plot(S_T[:, 0], S_T[:, 1], '-o')
    print N1, S_T[-1, 1]