# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:13:03 2020

@author: Alfre
"""

def interp_krg_cubique(x_obs, y_obs, z_obs, x_int, y_int, a, C):
    # Interpolation par splines
    # x_obs, y_obs, z_obs : observations
    # [np.array dimension 1*n]
    # x_int, y_int, positions pour lesquelles on souhaite interpoler une valeur z_int
    # [np array dimension m*p]
    # rho, pondération
    # int
    
    cote_A = x_obs.shape[0]
    z_int = np.nan*np.zeros(x_int.shape)
    A = np.ones((cote_A+1,cote_A+1))
    
    for i in np.arange(0,cote_A):
        for j in np.arange(0,cote_A):
            if i == j:
                A[i][j] = 0
                
            else:
                dist = np.sqrt((x_obs[i]-x_obs[j])**2+(y_obs[i]-y_obs[j])**2)
                A[i][j] = calculer_gamma0_cubique(C,a,dist)
    
    A[-1][-1] = 0
    
    # Construction des matrices A et B
    for k in np.arange(0,x_int.shape[0]):
        for l in np.arange(0,x_int.shape[1]):
            B = np.ones((cote_A+1,1))
            
            for i in np.arange(0,cote_A):
                dist_0 = np.sqrt((x_obs[i]-x_int[k,l])**2+(y_obs[i]-y_int[k,l])**2)
                B[i][0] = calculer_gamma0_cubique(C,a,dist_0)
            
            S = np.linalg.solve(A,B)
            z_int[k,l] = 0
            for i in np.arange(0,cote_A):
                z_int[k,l] += S[i]*z_obs[i][0] 
                
            
    return z_int