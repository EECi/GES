# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 09:18:07 2021

@author: rmw61
"""
import numpy as np
from dummy_parameters import T_k

def lamorturb(Gr, Re):
    
    Le = 0.819

    free = Gr < 1e5
    Nu_G = 0.5 * free * Gr**0.25 + 0.13*(1-free)*Gr**0.33

    forced = Re < 2e4
    Nu_R = 0.6*forced*Re**0.5 + 0.032*(1-forced)*Re**0.8  

    x = Nu_G > Nu_R

    Nu = x*Nu_G + (1-x)*Nu_R

    Sh = x*Nu*Le**0.25 + (1-x)*Nu*Le**0.33

    return(Nu, Sh)

def convection(d, A, T1, T2, ias):
    
    g = 9.81
    nu = 15.1e-6
    lam = 0.025
    
    Gr = (g*d**3)/(T1*nu**2)*abs(T1 - T2)
    Re = ias*d/nu
    (Nu, Sh) = lamorturb(Gr,Re)
    
    QV_1_2 = A*Nu*lam*(T1-T2)/d
    QP_1_2 = 0
    
    return(QV_1_2, QP_1_2)

def radiation(eps_1, eps_2, rho_1, rho_2, F_1_2, F_2_1, A_1, T_1, T_2):
    
    sigm = 5.67e-8
    
    k = eps_1*eps_2/(1-rho_1*rho_2*F_1_2*F_2_1)
    QR_1_2 = k*sigm*A_1*F_1_2*(T_1**4 - T_2**4)
    
    return(QR_1_2)

def conduction(A, lam, l, T1, T2):
    QD_12 = (A*lam/l)*(T1-T2)
    
    return(QD_12)
    
def T_ext(t):
    # Weather data

    climate = np.genfromtxt('climate.txt', delimiter=',')
    
    deltaT = 600
    n = int(np.ceil(t/deltaT))
    T_e = climate[n, 0] + T_k
    
    return(T_e)
    
def sat_conc(T):

    TC = T - T_k   
    spec_hum = np.exp(11.56 - 4030/(TC + 235))
    air_dens = -0.0046*TC + 1.2978
    a = spec_hum*air_dens
    
    return a
    
def Cw_ext(t):
    # Weather data

    climate = np.genfromtxt('climate.txt', delimiter=',')
    
    deltaT = 600
    n = int(np.ceil(t/deltaT))
    RH_e = climate[n, 1]/100;
    
    Cw_e = RH_e * sat_conc(T_ext(t))
    
    return(Cw_e)