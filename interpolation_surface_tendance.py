# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:54:02 2020

@author: MENGIN Alfred, GENET Arthur
"""

def interp_sfc(x_obs, y_obs, z_obs, x_int, y_int,p):
    # Interpolation par surfaces de tendance
    # x_obs, y_obs, z_obs : observations
    # [np.array dimension 1*n]
    # x_int, y_int, positions pour lesquelles on souhaite interpoler une valeur z_int
    # [np array dimension m*p]
    # p, pondération
    # int
    
    z_int = np.zeros(x_int.shape)
    
    nb_l_A = x_obs.shape[0]
    nb_c_A = int((p+2)*(p+1)/2)
    A = np.zeros((nb_l_A,nb_c_A))

    for k in np.arange(0,nb_l_A):
        compteur = 0
        for i in np.arange(0,p+1):
            for j in np.arange(0,p-i+1):
                A[k][compteur] = x_obs[k][0]**i*y_obs[k][0]**j
                compteur += 1
    alpha = np.dot(np.linalg.pinv(A),z_obs)
    
    for k in np.arange(0,x_int.shape[0]):
        for l in np.arange(0,x_int.shape[1]):
            compteur = 0
            for i in np.arange(0,p+1):
                for j in np.arange(0,p-i+1):
                    z_int[k,l] += alpha[compteur]*x_int[k,l]**i*y_int[k,l]**j
                    compteur += 1
    return z_int
