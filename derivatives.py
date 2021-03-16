# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 10:07:36 2021

@author: rmw61
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from functions import convection, radiation, conduction, T_ext, Cw_ext, sat_conc
from scipy.integrate import solve_ivp
from dummy_parameters import R, M_w, M_a, atm, H_fg, N_A, heat_phot
from dummy_parameters import V, A_c, A_f, A_v, A_m, A_p, A_l 
from dummy_parameters import d_c, d_f, d_m, d_p, cd_c, c_i, c_f, c_m, c_p
from dummy_parameters import F_c_f, F_f_c, F_c_v, F_c_m, F_l_c, F_l_v, F_l_m, F_l_p 
from dummy_parameters import F_m_l, F_f_p, F_c_l, F_m_v, F_v_l, F_p_l
from dummy_parameters import F_p_f, F_p_v, F_p_m, F_v_c, F_v_p, F_v_m, F_m_c, F_m_p 
from dummy_parameters import eps_c, eps_f, eps_v, eps_m, eps_p, eps_l
from dummy_parameters import rho_c, rho_f, rho_v, rho_m, rho_p, rho_l
from dummy_parameters import lam_c, l_c, rhod_c, c_c, lam_f, l_f, T_fl, lam_p, l_m
from dummy_parameters import T_ss, T_al
from dummy_parameters import ACH, ias, dehumidify
from dummy_parameters import f_heat, f_light, P_al, P_ambient_al
from dummy_parameters import c_v, msd_v, d_v, AF_g, LAI

tic = time.time()

# Initial conditions
T_i_0 = 300
T_c_0 = 295
T_f_0 = 295
T_v_0 = 295
T_m_0 = 295
T_p_0 = 295
T_c1_0 = 295
T_c2_0 = 290
T_c3_0 = 289
T_c4_0 = 288
T_c5_0 = 287
C_w_0 = 0.0054

z = [T_c_0, T_i_0, T_v_0, T_m_0, T_p_0, T_f_0, T_c1_0, T_c2_0, T_c3_0, T_c4_0, T_c5_0, C_w_0]
daynum = [0]

def day(t):
    ## Day
    day_new = np.ceil(t/86400)
    return(day_new)
    
def model(t,z):
    T_c = z[0]
    T_i = z[1]
    T_v = z[2]
    T_m = z[3]
    T_p = z[4]
    T_f = z[5]
    T_c1 = z[6]
    T_c2 = z[7]
    T_c3 = z[8]
    T_c4 = z[9]
    T_c5 = z[10]
    C_w = z[11]
    
    p_w = C_w*R*T_i/M_w
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i)
    
    daynum.append(day(t))
    if daynum[(len(daynum)-1)] > daynum[(len(daynum)-2)]:
        print(daynum[len(daynum)-1])
        
    ## Lights
    
    hour = np.floor(t/3600) + 1
    day_hour=(hour/24-np.floor(hour/24))*24
    L_on = (day_hour>-0.01 and day_hour<09.01) or day_hour > 15.01
    AL_on = day_hour>08.01 and day_hour<16.01
    
    T_l = L_on*T_al + (1-L_on)*T_i;
    
    QV_l_i = f_heat*P_al*L_on + P_ambient_al*AL_on
        
    ## Convection
    # Convection internal air -> cover

    (QV_i_c, QP_i_c) = convection(d_c, A_c, T_i, T_c, ias)

    # Convection internal air -> floor
    
    (QV_i_f, QP_i_f) = convection(d_f, A_f, T_i, T_f, ias)
    
    # Convection internal air -> vegetation
    A_v_exp = LAI*A_v
    (QV_i_v, QP_i_v) = convection(d_v, A_v_exp, T_i, T_v, ias)
    
    # Convection internal air -> mat
    A_m_exp = A_m*(1-AF_g)
    (QV_i_m, QP_i_m) = convection(d_m, A_m_exp, T_i, T_m, ias)
    
    # Convection internal air -> tray
    
    (QV_i_p, QP_i_p) = convection(d_p, A_p, T_i, T_p, ias)
    
    ## Radiation
    # Radiation cover to floor
    QR_c_f = radiation(eps_c, eps_f, rho_c, rho_f, F_c_f, F_f_c, A_c, T_c, T_f)
    
    # Radiation cover to vegetation
    QR_c_v = radiation(eps_c, eps_v, rho_c, rho_v, F_c_v, F_v_c, A_c, T_c, T_v)
    
    # Radiation cover to mat
    QR_c_m = radiation(eps_c, eps_m, rho_c, rho_m, F_c_m, F_m_c, A_c, T_c, T_m)
    
    # Radiation lights to cover
    QR_l_c = radiation(eps_l, eps_c, rho_l, rho_c, F_l_c, F_c_l, A_l, T_l, T_c)
    
    # Radiation lights to vegetation
    QR_l_v = radiation(eps_l, eps_v, rho_l, rho_v, F_l_v, F_v_l, A_l, T_l, T_v)
    
    # Radiation lights to mat
    QR_l_m = radiation(eps_l, eps_m, rho_l, rho_m, F_l_m, F_m_l, A_l, T_l, T_m)
    
    # Radiation lights to tray
    QR_l_p = radiation(eps_l, eps_p, rho_l, rho_p, F_l_p, F_p_l, A_l, T_l, T_p)
    
    # Radiation vegetation to cover
    QR_v_c = radiation(eps_v, eps_c, rho_v, rho_c, F_v_c, F_c_v, A_v, T_v, T_c)
    
    # Radiation vegetation to mat
    QR_v_m = radiation(eps_v, eps_m, rho_v, rho_m, F_v_m, F_m_v, A_v, T_v, T_m)
    
    # Radiation vegetation to tray
    QR_v_p = radiation(eps_v, eps_p, rho_v, rho_p, F_v_p, F_p_v, A_v, T_v, T_p)
    
    # Radiation mat to cover
    QR_m_c = radiation(eps_m, eps_c, rho_m, rho_c, F_m_c, F_c_m, A_m, T_m, T_c)
    
    # Radiation mat to vegetation
    QR_m_v = radiation(eps_m, eps_v, rho_m, rho_v, F_m_v, F_v_m, A_m, T_m, T_v)
    
    # Radiation mat to tray
    QR_m_p = radiation(eps_m, eps_p, rho_m, rho_p, F_m_p, F_p_m, A_m, T_m, T_p)
    
    # Radiation tray to vegetation
    QR_p_v = radiation(eps_p, eps_v, rho_p, rho_v, F_p_v, F_v_p, A_p, T_p, T_v)
    
    # Radiation tray to mat
    QR_p_m = radiation(eps_p, eps_m, rho_p, rho_m, F_p_m, F_m_p, A_p, T_p, T_m)
    
    # Radiation tray to floor
    QR_p_f = radiation(eps_p, eps_f, rho_p, rho_f, F_p_f, F_f_p, A_p, T_p, T_f)
       
    # Radiation floor to cover
    QR_f_c = radiation(eps_f, eps_c, rho_f, rho_c, F_f_c, F_c_f, A_f, T_f, T_c)
    
    # Radiation floor to tray
    QR_f_p = radiation(eps_f, eps_p, rho_f, rho_p, F_f_p, F_p_f, A_f, T_f, T_p)
    
    
    ## Conduction
    # Conduction through cover
    QD_c12 = conduction(A_c, lam_c[0], l_c[0], T_c, T_c1)
    QD_c23 = conduction(A_c, lam_c[1], l_c[1], T_c1, T_c2)
    QD_c34 = conduction(A_c, lam_c[2], l_c[2], T_c2, T_c3)
    QD_c45 = conduction(A_c, lam_c[3], l_c[3], T_c3, T_c4)
    QD_c56 = conduction(A_c, lam_c[4], l_c[4], T_c4, T_c5)
    QD_c67 = conduction(A_c, lam_c[5], l_c[5], T_c5, T_ss)
    
    # Conduction through floor
    QD_f12 = conduction(A_f, lam_f, l_f, T_f, T_fl)
    
    # Conduction mat to tray
    QD_m_p = (A_m*lam_p/l_m)*(T_m-T_p)
    
    ## Transpiration
    QS_int = f_light*P_al*L_on/A_p

    PPFD = QS_int/1e-6/N_A/heat_phot
    r_aG = 100
    r_sG = 60*(1500+PPFD)/(200+PPFD)
    QT_G = A_v*(2*LAI*H_fg*(1/(r_aG+r_sG))*(sat_conc(T_v) - C_w));

    QT_v_i = QT_G
    
    ## Ventilation
    
    QV_i_e = ACH*V*rho_i*c_i*(T_i - T_ext(t))
    
    MW_i_e = ACH*(C_w - Cw_ext(t))
    
    ## Dehumidification

    MW_cc_i = -1*dehumidify/3600
    
    #print(MW_i_e)
    
    # ODE equations
    
    dT_cdt = (1/(A_c*cd_c))*(QV_i_c - QR_c_f - QR_c_v - QR_c_m
                             + QR_l_c - QD_c12)
    dT_idt = (1/(V*rho_i*c_i))*(- QV_i_c - QV_i_f - QV_i_e + QV_l_i - QV_i_m 
                                - QV_i_v - QV_i_p)
    dT_fdt = (1/(A_f*c_f))*(QV_i_f - QR_f_c - QR_f_p - QD_f12)
    dT_vdt = (1/(c_v*A_v*msd_v))*(QV_i_v - QR_v_c - QR_v_m + QR_l_v - QR_v_p - QT_v_i)
    dT_mdt = (1/(A_m*c_m))*(QV_i_m + QP_i_m - QR_m_v - QR_m_c + QR_l_m - QR_m_p - QD_m_p)
    dT_pdt = (1/(A_p*c_p))*(QD_m_p + QV_i_p + QP_i_p - QR_p_f + QR_l_p - QR_p_v - QR_p_m)
    dT_c1dt = (1/(rhod_c[1]*c_c[1]*l_c[1]*A_c))*(QD_c12-QD_c23)
    dT_c2dt = (1/(rhod_c[2]*c_c[2]*l_c[2]*A_c))*(QD_c23-QD_c34)
    dT_c3dt = (1/(rhod_c[3]*c_c[3]*l_c[3]*A_c))*(QD_c34-QD_c45)
    dT_c4dt = (1/(rhod_c[4]*c_c[4]*l_c[4]*A_c))*(QD_c45-QD_c56)
    dT_c5dt = (1/(rhod_c[5]*c_c[5]*l_c[5]*A_c))*(QD_c56-QD_c67)
    
    dC_wdt = (1/(V*H_fg))*(QT_v_i-QP_i_c-QP_i_f-QP_i_m-QP_i_p) - MW_i_e + MW_cc_i
    #dC_wdt = -MW_i_e
    
    return np.array([dT_cdt,dT_idt,dT_vdt,dT_mdt,dT_pdt,dT_fdt,dT_c1dt,
                     dT_c2dt,dT_c3dt,dT_c4dt,dT_c5dt,dC_wdt])


t = [0,432000]
tval = np.linspace(0,432000,121)

output = solve_ivp(model, t, z, method='BDF', t_eval=tval, rtol = 1e-3)

toc = time.time()

print(toc-tic)

# Plot Outputs

Temperatures = np.transpose(output.y[:11,:])
MoistureContent = np.transpose(output.y[11,:])
plt.figure(1)
plt.plot(output.t,Temperatures)
plt.figure(2)
plt.plot(output.t,MoistureContent)