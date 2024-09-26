"""

#===============================================================================
#                 INTERPOLATE A SURFACE Z(X,Y) WHERE MAP(X,Y)=1                #
#===============================================================================

Interpolate a boolean_map(x,y) to boolean_map_new(xnew,ynew) with xnew = step*x 
and ynew = step*y. The boolean_map is generally used to piece-wise interpolate 
z(x,y) surfaces, where boolean_map(x, y) states which points are kept or not.

If <removeBoundaries> = True, we make sure that we don't extrapolate towards
points of which we have no information. If we have a map z(x,y) where we have 
a correlated piece left of x=5 and one right of x=5, then at x=5 we should have
a black boundary since we don't know what happens between x=[4,6] for example.

Hanne Thienpondt
01/09/2022

"""

#!/usr/bin/python3 
import numpy as np 

#===============================================================================
#                 INTERPOLATE A SURFACE Z(X,Y) WHERE MAP(X,Y)=1                #
#===============================================================================

def interpolate_booleanMap(x,y,xnew,ynew,boolean_map,removeBoundaries=True):  
    
    # Return immediately if the boolean_map only contains True's or False's
    if np.all(boolean_map==False): return np.zeros((len(xnew),len(ynew))).astype(np.bool) 
    if np.all(boolean_map==True):  return np.ones((len(xnew),len(ynew))).astype(np.bool)  
    
    # Intiate <boolean_map_new>(xnew,ynew) at a higher resolution than <boolean_map>(x,y)
    boolean_map_new = np.ones((len(xnew),len(ynew))) 
           
    # If the original square was removed, then remove it again: boolean_map(x1,y1) = 0 --> boolean_map_new(x2:x3,y2:y3) = 0
    # If we increase from x=[0,1,2] to xnew=[0,0.2,0.4,0.6,0.8,1,...2] then tile [1] becomes [0.6,0.8,1.0,1.2,1.4]
    # So xnew ranges from xold_bot=0.5 to xold_top=1.5 if x=1
    xold_bot = [x[0]] + [x[i] - (x[i]-x[i-1])/2 for i in range(len(x)) if i!=0]
    xold_top = [x[i] + (x[i+1]-x[i])/2 for i in range(len(x)-1)] + [x[-1]]
    yold_bot = [y[0]] + [y[i] - (y[i]-y[i-1])/2 for i in range(len(y)) if i!=0]
    yold_top = [y[i] + (y[i+1]-y[i])/2 for i in range(len(y)-1)] + [y[-1]] 
    for i in range(len(x)):   
        for j in range(len(y)):   
            indexes_k = [k for k in range(len(xnew)) if xnew[k]>=xold_bot[i] and xnew[k]<=xold_top[i]]
            indexes_l = [k for k in range(len(ynew)) if ynew[k]>=yold_bot[j] and ynew[k]<=yold_top[j]]  
            boolean_map_new[np.ix_(indexes_k,indexes_l)] = boolean_map[i,j]  
            
    # Return the new map if we don't care about the boundaries between True's and False's 
    if not removeBoundaries: return boolean_map_new.astype(np.bool)
    
    # The <boolean_map>(x, y) states which points are retained or not
    # If we have for x=0 and y=[0,1,2,3]: True, True, False, False
    # Then the boundary between the two areas lies at (x=0, y=1.5)
    boundaries = []
    
    # Collect the indices of the boundaries 
    for row in range(np.shape(boolean_map)[0]):
        if True in boolean_map[row,:] and False in boolean_map[row,:]:
            for i in range(len(boolean_map[row,:])-1):
                if boolean_map[row,i]!=boolean_map[row,i+1]: boundaries.append([row, i+0.5]) 
    for column in range(np.shape(boolean_map)[1]):
        if True in boolean_map[:,column] and False in boolean_map[:,column]:
            for i in range(len(boolean_map[:,column])-1):
                if boolean_map[i,column]!=boolean_map[i+1,column]: boundaries.append([i+0.5,column]) 

    # Iterate through each square of the boolean_map, e.g. x=[3,4] and y=[0,1] is a square
    # If (x=3.5, y=0) and (x=3.5, y=1) are boundaries, remove the entire square
    # If (x=3, y=0.5) and (x=4, y=0.5) are boundaries, remove the entire square
    for i in range(len(x)-1): 
        for j in range(len(y)-1):  
               
            # Check the boundaries for the square [i:i+1, j,j+1]
            if (([i, j+0.5] in boundaries) and ([i+1, j+0.5] in boundaries))\
            or (([i+0.5, j] in boundaries) and ([i+0.5, j+1] in boundaries)):
               
                # The square in <boolean_map> ranges from [i:i+1, j,j+1]
                # Find the indices [kstart:kend, lstart:lend] of the square in <boolean_map_new>
                indexes_k = [k for k in range(len(xnew)) if xnew[k]>=x[i] and xnew[k]<x[i+1]]
                indexes_l = [k for k in range(len(ynew)) if ynew[k]>=y[j] and ynew[k]<y[j+1]] 
                   
                # If the square needs to be removed, indicate it with <new_boolean_map>[k,l] = 1
                for k in indexes_k:
                    for l in indexes_l:
                        boolean_map_new[k,l] = 0
       
    # If (x=3.5, y=0) and (x=4, y=0.5) are boundaries, remove half the square (triangle) 
    for i in range(len(x)-1): 
        for j in range(len(y)-1):  
                
            # Check the boundaries for the square [i:i+1, j,j+1]
            if (([i+0.5, j] in boundaries)   and ([i+1, j+0.5] in boundaries))\
            or (([i+1, j+0.5] in boundaries) and ([i+0.5, j+1] in boundaries))\
            or (([i+0.5, j+1] in boundaries) and ([i, j+0.5] in boundaries))\
            or (([i, j+0.5] in boundaries)   and ([i+0.5, j] in boundaries)):
                
                # The square in <boolean_map> ranges from [i:i+1, j,j+1]
                # Find the indices [kstart:kend, lstart:lend] of the square in <boolean_map_new>
                indexes_k = [k for k in range(len(xnew)) if xnew[k]>=x[i] and xnew[k]<x[i+1]]
                indexes_l = [l for l in range(len(ynew)) if ynew[l]>=y[j] and ynew[l]<y[j+1]]  
                    
                # If the half of the square needs to be removed, indicate it with <new_boolean_map>[k,l] = 1
                for k in indexes_k:
                    for l in indexes_l:
                        xi = round(xnew[k]-x[i], 10)
                        yi = round(ynew[l]-y[j], 10)
                        if ([i+0.5, j] in boundaries) and ([i+1, j+0.5] in boundaries): 
                            if yi <= xi: boolean_map_new[k,l] = 0
                        if ([i, j+0.5] in boundaries) and ([i+0.5, j] in boundaries):  
                            if yi+xi < 1: boolean_map_new[k,l] = 0
                        if ([i+1, j+0.5] in boundaries) and ([i+0.5, j+1] in boundaries): 
                            if yi+xi >= 1: boolean_map_new[k,l] = 0
                        if ([i+0.5, j+1] in boundaries) and ([i, j+0.5] in boundaries): 
                            if yi >= xi: boolean_map_new[k,l] = 0
  
    return boolean_map_new.astype(np.bool)
