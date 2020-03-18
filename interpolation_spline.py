# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:52:18 2020

@author: MENGIN Alfred, GENET Arthur
"""

def phi(x_obs, y_obs, x_int, y_int):
    dist = np.sqrt((x_int-x_obs)**2+(y_int-y_obs)**2)
    return dist**2*np.log(dist)


def interp_spl(x_obs, y_obs, z_obs, x_int, y_int, rho):
    # Interpolation par splines
    # x_obs, y_obs, z_obs : observations
    # [np.array dimension 1*n]
    # x_int, y_int, positions pour lesquelles on souhaite interpoler une valeur z_int
    # [np array dimension m*p]
    # rho, pondération
    # int
    
    cote_A = x_obs.shape[0]
    z_int = np.nan*np.zeros(x_int.shape)
    A = np.zeros((cote_A+3,cote_A+3))
    B = np.zeros((cote_A+3,1))
    
    # Construction des matrices A et B
    for i in np.arange(0,cote_A):
        A[i][0] = 1
        A[i][1] = x_obs[i][0]
        A[i][2] = y_obs[i][0]
        A[cote_A][i+3] = 1
        A[cote_A+1][i+3] = x_obs[i][0]
        A[cote_A+2][i+3] = y_obs[i][0]
        for j in np.arange(0,cote_A):
            if i == j:
                A[i][j+3] = rho
                
            else:
                A[i][j+3] = phi(x_obs[i][0],y_obs[i][0],x_obs[j][0],y_obs[j][0])
                
        B[i][0] = z_obs[i][0]        
    S = np.linalg.solve(A,B)
    
    for i in np.arange(0,x_int.shape[0]):
        for j in np.arange(0,x_int.shape[1]):
            z_int[i,j] = S[0]  + S[1]* x_int[i,j] + S[2] * y_int[i,j]
            for k in np.arange(0,cote_A):
                z_int[i,j] += S[k+3]*phi(x_obs[k][0],y_obs[k][0],x_int[i,j],y_int[i,j])
    return z_int