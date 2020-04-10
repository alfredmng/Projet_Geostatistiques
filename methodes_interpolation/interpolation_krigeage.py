# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:13:03 2020

@author: Alfre
"""

import numpy as np

def calculer_gamma0_cubique(C,a,h):
    return C*((7*h**2)/a**2 - (8.75*h**3)/a**3 + (3.5*h**5)/a**5 - (0.75*h**7)/a**7)

def calculer_deriveC0_cubique(C,a,h):
    return (7*h**2)/a**2 - (8.75*h**3)/a**3 + (3.5*h**5)/a**5 - (0.75*h**7)/a**7

def calculer_derivea0_cubique(C,a,h):
    return C*(-(14*h**2)/a**3 + (26.25*h**3)/a**4 - (17.5*h**5)/a**6 + (5.25*h**7)/a**8)
            
def construire_A_B_cubique(nuee,C,a):
    
    A = []
    B = []
    for i in np.arange(0,nuee.shape[0]):
        
        h = nuee[i][3]
        
        ligne_A = []
        ligne_A.append(calculer_derivea0_cubique(C,a,h))
        ligne_A.append(calculer_deriveC0_cubique(C,a,h))
        A.append(ligne_A)
        
        gamma0 = calculer_gamma0_cubique(C,a,h)
        gammaexp = nuee[i][2]
        ligne_B = []
        ligne_B.append(gammaexp-gamma0)
        B.append(ligne_B)
        
    return [np.asarray(A),np.asarray(B)]



def estimation_a_C_cubique(nuee,ecarta,ecartC,rang_max):
    
    Csuiv = trouver_a0_C0(nuee)[0]
    asuiv = trouver_a0_C0(nuee)[1]
    
    
    Cprec = 10**8
    aprec = 10**8
    
    rang= 0
    res = []
    
    while abs(Csuiv-Cprec) > ecartC and abs(asuiv - aprec) > ecartC and rang < rang_max:
      
        rang += 1
        
        A = construire_A_B_cubique(nuee,Csuiv,asuiv)[0]
        
        B = construire_A_B_cubique(nuee,Csuiv,asuiv)[1]
        P = np.eye(A.shape[0]) 
        
        X_chap = (moindres_carres(A,B,P)[0]).reshape(2,)
        V_chap = moindres_carres(A,B,P)[1]
        
        
        # Mise à jour
        garage = asuiv
        aprec = asuiv
        asuiv = garage + X_chap[0]
        
        garage = Csuiv
        Cprec = Csuiv
        Csuiv = garage + X_chap[1]
        
        ligne = [asuiv,Csuiv]
        res.append(ligne)
        
        print("Itération " + str(rang))
    return res[-1]

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