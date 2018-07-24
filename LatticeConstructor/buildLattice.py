import numpy as np
import pandas as pd

def getXYCoords(n1,n2, a1_x, a1_y, a2_x, a2_y):
    x = n1*a1_x + n2*a2_x
    y = n1*a1_y + n2*a2_y
    return (x,y)
    


def buildLattice(N1, N2, alpha, neigb_cutOff, 
                 a1 = 1, 
                 a2 = 1, 
                 w0 = 1, 
                 power = 3):
    precision = 3
    # the lattice primitive vectors
    # a1 always  = x
    a1_x = 1 * a1
    a1_y = 0 * a1
    #
    a2_x = np.round(np.cos(alpha) * a2, precision)
    a2_y = np.round(np.sin(alpha) * a2, precision)
    #
#    print a1_x, a1_y
#    print a2_x, a2_y
    
    neigbs_table = []
    weight_table = []
    for i2 in range(N2):
        for i1 in range(N1):
            neigbs_column = []
            neigbs_column.append(i2*N1 + i1)
            
            weight_column = []
            weight_column.append(0)
            
            # coords of i1, i2
            (x0, y0) = getXYCoords(i1,i2, a1_x, a1_y, a2_x, a2_y)
            
            for j2 in range(0, N2):
                for j1 in range(0, N1):
                    
                    if j1==i1 and j2==i2:
                        pass
                    else:
                        (x_n, y_n) = getXYCoords(j1, j2, a1_x, a1_y, a2_x, a2_y)
                        distance = np.sqrt((x0 - x_n)**2 + (y0 - y_n)**2)
                        ### periodic a1
                        # cartesian coords
                        x=x0+N1*a1_x
                        y=y0+N1*a1_y
                        distance10 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        #
                        x=x0-N1*a1_x
                        y=y0-N1*a1_y
                        distance11 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        
                        distance1=np.min([distance10, distance11])

                        ### peridic in a2
                        x=x0+N2*a2_x
                        y=y0+N2*a2_y
                        distance20 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        #
                        x=x0-N2*a2_x
                        y=y0-N2*a2_y
                        distance21 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        distance2=np.min([distance20, distance21])
                        ### periodic in both a1 and a2:
                        ### peridic in a2
                        x=x0 + N1*a1_x + N2*a2_x
                        y=y0 + N1*a1_y + N2*a2_y
                        distance30 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        #
                        x=x0 + N1*a1_x - N2*a2_x
                        y=y0 + N1*a1_y - N2*a2_y
                        distance31 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        #
                        x=x0 - N1*a1_x + N2*a2_x
                        y=y0 - N1*a1_y + N2*a2_y
                        distance32 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        #
                        x=x0 - N1*a1_x - N2*a2_x
                        y=y0 - N1*a1_y - N2*a2_y
                        distance33 = np.sqrt((x - x_n)**2 + (y - y_n)**2)
                        distance3=np.min([distance30, 
                                          distance31,
                                          distance32,
                                          distance33,
                                          ])                        

                        

                        distance = np.min([
                                       distance, 
                                       distance1,
                                       distance2,
                                       distance3,
                                       ])
                        ####
                        neigbs_column.append(j2*N1 + j1)
                        # weight is 1./|distance|**3 * w0
                        # we can randomly replace F <-> AF <==> w0 <--> -W0
                        weight = np.round(w0 * 1/((distance +0.0)**power), precision)
                        weight_column.append(weight)
            ###### save neigbs/weights  of  i1, i2
            weight_table.append(weight_column)                
            neigbs_table.append(neigbs_column)
    # make a table            
    weight_array = np.array(weight_table)                
    neighbs_array = np.array(neigbs_table)                
    ################################################
    w_df = pd.DataFrame(weight_array)
    w_df = w_df.transpose()
    n_df = pd.DataFrame(neighbs_array)
    n_df = n_df.transpose()
    full_df = pd.concat([w_df, n_df], axis=1)
#    print full_df
#    # get the neigbors unique weight up to the cutoff
    unique_neigbs_upto_cutOff = np.array(list(set(full_df[0].values[:, 0])))
    unique_neigbs_upto_cutOff.sort()
    #
    if neigb_cutOff>len(unique_neigbs_upto_cutOff):
        neigb_cutOff = len(unique_neigbs_upto_cutOff)
    #
    unique_neigbs_upto_cutOff = unique_neigbs_upto_cutOff[::-1][:neigb_cutOff]
     # to include the site i1, i2
    unique_neigbs_upto_cutOff = np.array([0] + list(unique_neigbs_upto_cutOff))
#    print unique_neigbs_upto_cutOff
#    # extract the relevant part
    weights_dff = pd.DataFrame()
    neighbs_dff = pd.DataFrame()
    for i in range(N1*N2):
        dff = pd.DataFrame()
        for n in range(neigb_cutOff+1):
            weight = unique_neigbs_upto_cutOff[n]
            df = full_df.loc[full_df[i].values[:, 0] == weight]
    #        print df 
            dff = pd.concat([dff, df], axis=0)
            
        weights = dff.values[:, i]
        weights_df = pd.DataFrame(weights)
    #    print weights_df.shape
        weights_dff = pd.concat([weights_dff, weights_df], axis=1)
#        print weights_dff
        
        neighbs = dff.values[:, N1*N2 + i]  
        neighbs_df = pd.DataFrame(neighbs)
    #    print weights_df.shape
        neighbs_dff = pd.concat([neighbs_dff, neighbs_df], axis=1)
#        print neighbs_dff
#    print neighbs
    return weights_dff.values[:, :], neighbs_dff.values[:,:]
    
    
######################################################################

    
##### test if all the neighbs are correct    
if __name__=='__main__':
    import pandas as pd
    N1 = 4
    N2 = 4
    neigb_cutOff = 3
    theta_factor = 3 #
    alpha = np.pi/ theta_factor
    w_clean, n_clean = buildLattice(N1, N2, alpha, neigb_cutOff)
##    w_unclean, n_unclean, w, n = buildLattice(N1, N2, alpha, neigb_cutOff)    
#    w_df = pd.DataFrame(w)
#    w_df = w_df.transpose()
##    w_df.loc[w_df[i][:] = ]
#    unique_neigbs_upto_neigb_cutOff = []
#    
#    for i in range(w_df.shape[1]):
#        q = np.array(list(set(w_df[i][:])))
#        q.sort()
#        
#        unique_neigbs_upto_neigb_cutOff.append(q[::-1][:neigb_cutOff])
#    
##    unique_neighbs = unique_neigbs_upto_neigb_cutOff[6][:5]
##    for neighb in unique_neighbs:
##        print ' **********************************   ' , neighb
##        for column in range(N1*N2):
##            print column, ' --- > ' , w_df.loc[ w_df[column][:]==neighb ].shape[0]
#    
#    ###how to attach the frames and clean for up_to_rank
#    ##################################################################
#    n_df = pd.DataFrame(n)
#    n_df = n_df.transpose()
#    full_df = pd.concat([w_df, n_df], axis=1)
#    # get the neigbors unique weight up to the cutoff
#    unique_neigbs_upto_cutOff = np.array(list(set(full_df[0].values[:, 0])))
#    unique_neigbs_upto_cutOff.sort()
#    unique_neigbs_upto_cutOff = unique_neigbs_upto_cutOff[::-1][:neigb_cutOff]
#    # to include the site i1, i2
#    unique_neigbs_upto_cutOff = np.array([0] + list(unique_neigbs_upto_cutOff))
#    # extract the relevant part
#    dff = pd.DataFrame()
#    for n in range(neigb_cutOff+1):
#        weight = unique_neigbs_upto_cutOff[n]
#        df = full_df.loc[full_df[0].values[:, 0] == weight]
##        print df 
#        dff = pd.concat([dff, df], axis=0)
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
    