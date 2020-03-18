# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:08:47 2020

@author: MENGIN Alfred, GENET Arthur
"""

####################### Fonctions auxiliaires ######################

def trouver_a0_C0(nuee):
    C0 = 0
    a0 = 0
    for i in np.arange(0,nuee.shape[0]):
        if nuee[i][2] > C0:
            C0 = nuee[i][2]
            a0 = nuee[i][3]
    return [C0,a0]

def moindres_carres(A,B,P):
    
    """ 
    Calcule la matrice des paramètres et le vecteur des résidus
    
    :param A : la matrice modèle
    :type A : array
    :param B : la matrice des pseudo_observations
    :type B : array
    :param P : matrice de poids
    :type P : array
    
    :return : la liste contenant le vecteur des paramètres estimé et le vecteur des résidus 
    :rtype : list(array) 
    """
    
    K = np.dot(np.dot(A.transpose(),P),B)
    N = np.dot(np.dot(A.transpose(),P),A)
    if np.linalg.det(N) != 0:
        X_chap = np.dot(np.linalg.inv(N),K)
        V_chap = B - np.dot(A,X_chap)
        return [X_chap,V_chap]
    
####################### Variogramme analytique cubique ######################
    
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
        
def variogramme_ana_cubique(nuee,ecarta,ecartC,rang_max):
    abscisses = np.linspace(0,200,200)
    gamma_ana = []
    param = estimation_a_C_cubique(nuee,ecarta,ecartC,rang_max)
    a = param[0]
    C = param[1]
    for i in np.arange(0,abscisses.shape[0]):
        h = abscisses[i]
        if h < a:
            gamma_ana.append(calculer_gamma0_cubique(C,a,h))
        else:
            gamma_ana.append(C)
    
    plt.plot(abscisses,gamma_ana)
    #plt.scatter(tab_dist.transpose()[3],tab_dist.transpose()[2])
    plt.show()
    return gamma_ana  