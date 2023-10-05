# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 09:18:07 2021

@author: rmw61
"""
import numpy as np
from parameters import *

def lamorturb(Gr, Re):

    free = Gr < 1e5
    Nu_G = 0.5 * free * Gr**0.25 + 0.13*(1-free)*Gr**0.33

    forced = Re < 2e4
    Nu_R = 0.6*forced*Re**0.5 + 0.032*(1-forced)*Re**0.8  

    x = Nu_G > Nu_R

    Nu = x*Nu_G + (1-x)*Nu_R

    Sh = x*Nu*Le**0.25 + (1-x)*Nu*Le**0.33

    return(Nu, Sh)

def convection(d, A, T1, T2, ias, rho, c, C):
    
    g = 9.81
    nu = 15.1e-6
    lam = 0.025
    
    Gr = (g*d**3)/(T1*nu**2)*abs(T1 - T2)
    Re = ias*d/nu
    (Nu, Sh) = lamorturb(Gr,Re)
    
    QV_1_2 = A*Nu*lam*(T1-T2)/d
    QP_1_2 = A*H_fg/(rho*c)*Sh/Le*lam/d*(C - sat_conc(T2))
    #QP_1_2 = 0
    
    return(QV_1_2, QP_1_2, Nu)

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

def day(t):
    ## Day
    day_new = np.ceil(t/86400)
    return(day_new)
    
def model(t,z, climate, daynum):

    # Values being calculated

    T_c = z[0]
    T_i = z[1]
    T_v = z[2]
    T_m = z[3]
    T_p = z[4]
    T_f = z[5]
    T_s1 = z[6]
    T_s2 = z[7]
    T_s3 = z[8]
    T_s4 = z[9]
    T_vmean = z[10]
    T_vsum = z[11]
    C_w = z[12]
    C_c = z[13]
    C_buf = z[14]
    C_fruit = z[15]
    C_leaf = z[16]
    C_stem = z[17]
    R_fruit = z[18]
    R_leaf = z[19]
    R_stem = z[20]

    # External weather and dependent internal parameter values
    n = int(np.ceil(t/deltaT)) # count
    T_ext = climate[n, 0] + T_k # External air temperature (K)
    T_sk = climate[n, 1] + T_k # External sky temperature (K)
    wind_speed = climate[n, 2] # External wind speed (m/s)
    RH_e = climate[n, 3]/100 # External relative humidity
    Cw_ext = RH_e * sat_conc(T_ext) # External air moisture content
    p_w = C_w*R*T_i/M_w # Partial pressure of water [Pa]
    rho_i = ((atm - p_w)*M_a + p_w*M_w)/(R*T_i) # Internal density of air [kg/m^3]   
    LAI = SLA*C_leaf # Leaf area index
    C_ce = 4.0e-4*M_c*atm/(R*T_ext) # External carbon dioxide concentration [kg/m^3]
    C_c_ppm = C_c*R*T_i/(M_c*atm)*1.e6 # External carbon dioxide concentration [ppm]

    daynum.append(day(t)) # Day number

    # Option for printing progress of run in days - uncomment out if needed
    if daynum[(len(daynum)-1)] > daynum[(len(daynum)-2)]:
        print('Day', daynum[len(daynum)-1])
    
    hour = np.floor(t/3600) + 1
    # Option for printing progress in hours - uncomment if needed
    #print('Hour', hour)

    day_hour=(hour/24-np.floor(hour/24))*24

    ## Lights
    L_on = 0 # No additional lighting included
    AL_on = 0 # No ambient lighting included
        
    ## Convection
    # Convection external air -> cover

    (QV_e_c, QP_e_c, Nu_e_c ) = convection(d_c, A_c, T_ext, T_c, wind_speed, rho_i, c_i, C_w)
    QP_e_c = 0 # Assumed no external condensation/evaporation

    # Convection internal air -> cover

    (QV_i_c, QP_i_c, Nu_i_c) = convection(d_c, A_c, T_i, T_c, ias, rho_i, c_i, C_w)
    QP_i_c = max(QP_i_c,0) # assumed no evaporation from the cover, only condensation

    # Convection internal air -> floor
    
    (QV_i_f, QP_i_f, Nu_i_f) = convection(d_f, A_f, T_i, T_f, ias, rho_i, c_i, C_w)
    QP_i_f = max(QP_i_f,0) # assumed no evaporation from the floor, only condensation
    
    # Convection internal air -> vegetation
    A_v_exp = LAI*A_v
    (QV_i_v, QP_i_v, Nu_i_v) = convection(d_v, A_v_exp, T_i, T_v, ias, rho_i, c_i, C_w)
    QP_i_v = 0 # No condensation/evaporation - transpiration considered separately
    HV = Nu_i_v*lam/d_v
    
    # Convection internal air -> mat
    A_m_wool = 0.75*A_m # Area of mat exposed
    A_m_water = 0.25*A_m # assumed 25% saturated

    (QV_i_m, QP_i_m, Nu_i_m) = convection(d_m, A_m_wool, T_i, T_m, ias, rho_i, c_i, C_w)
    
    QP_i_m = A_m_water/A_m_wool * QP_i_m # Factored down
    
    # Convection internal air -> tray
    (QV_i_p, QP_i_p, Nu_i_p) = convection(d_p, A_p, T_i, T_p, ias, rho_i, c_i, C_w)
    QP_i_p = 0 # Assumed no condensation/evaporation from tray

    ## Far-IR Radiation

    A_vvf=min(LAI*p_v*A_f,p_v*A_f)
    F_c_v = min((1-F_c_f)*LAI,(1-F_c_f)) # Cover to vegetation
    F_c_m = max((1-F_c_f)*(1-LAI),0) # Cover to mat
    F_m_c = max((1-LAI),0.0) # Mat to cover
    F_m_v = 1-F_m_c # Mat to vegetation

    # Cover to sky
    QR_c_sk = radiation(eps_ce, 1, 0, 0, 1, 0, A_c, T_c, T_sk)

    # Radiation cover to floor
    QR_c_f = radiation(eps_ci, eps_s, rho_ci, rho_s, F_c_f, F_f_c, A_c_roof, T_c, T_f)
    
    # Radiation cover to vegetation
    QR_c_v = radiation(eps_ci, eps_v, rho_ci, rho_v, F_c_v, F_v_c, A_c_roof, T_c, T_v)
    
    # Radiation cover to mat
    QR_c_m = radiation(eps_ci, eps_m, rho_ci, rho_m, F_c_m, F_m_c, A_c_roof, T_c, T_m)
    
    # Radiation vegetation to cover
    QR_v_c = radiation(eps_v, eps_ci, rho_v, rho_ci, F_v_c, F_c_v, A_vvf, T_v, T_c)
    
    # Radiation vegetation to mat
    QR_v_m = radiation(eps_v, eps_m, rho_v, rho_m, F_v_m, F_m_v, A_vvf, T_v, T_m)
    
    # Radiation vegetation to tray
    QR_v_p = radiation(eps_v, eps_p, rho_v, rho_p, F_v_p, F_p_v, A_vvf, T_v, T_p)
    
    # Radiation mat to cover
    QR_m_c = radiation(eps_m, eps_ci, rho_m, rho_ci, F_m_c, F_c_m, A_m, T_m, T_c)
    
    # Radiation mat to vegetation
    QR_m_v = radiation(eps_m, eps_v, rho_m, rho_v, F_m_v, F_v_m, A_m, T_m, T_v)
    
    # Radiation mat to tray
    QR_m_p = radiation(eps_m, eps_p, rho_m, rho_p, F_m_p, F_p_m, A_m, T_m, T_p)
    
    # Radiation tray to vegetation
    QR_p_v = radiation(eps_p, eps_v, rho_p, rho_v, F_p_v, F_v_p, A_p, T_p, T_v)
    
    # Radiation tray to mat
    QR_p_m = radiation(eps_p, eps_m, rho_p, rho_m, F_p_m, F_m_p, A_p, T_p, T_m)
    
    # Radiation tray to floor
    QR_p_f = radiation(eps_p, eps_s, rho_p, rho_s, F_p_f, F_f_p, A_p, T_p, T_f)
       
    # Radiation floor to cover
    QR_f_c = radiation(eps_s, eps_ci, rho_s, rho_ci, F_f_c, F_c_f, A_f, T_f, T_c)
    
    # Radiation floor to tray
    QR_f_p = radiation(eps_s, eps_p, rho_s, rho_p, F_f_p, F_p_f, A_f, T_f, T_p)
    
    
    ## Conduction
    # Conduction through floor
    QD_sf1 = conduction(A_f, lam_s[0], l_s[0], T_f, T_s1)
    QD_s12 = conduction(A_f, lam_s[1], l_s[1], T_s1, T_s2)
    QD_s23 = conduction(A_f, lam_s[2], l_s[2], T_s2, T_s3)
    QD_s34 = conduction(A_f, lam_s[3], l_s[3], T_s3, T_s4)
    QD_s45 = conduction(A_f, lam_s[4], l_s[4], T_s4, T_ss)
    
    # Conduction mat to tray
    QD_m_p = (A_m*lam_p/l_m)*(T_m-T_p)

    ## Ventilation
    # Leakage (equations for orifice flow from Awbi, Ventilation of Buildings, Chapter 3)
    wind_speed_H = wind_speed*c*H**a # Wind speed at height H
    wind_pressure = Cp*0.5*rho_i*wind_speed_H**2 # Equals DeltaP for wind pressure
    stack_pressure_diff = rho_i*g*H*(T_i - T_ext)/T_i # DeltaP for stack pressure
 
    Qw = Cd*crack_area*(2*wind_pressure/rho_i)**0.5 # Flow rate due to wind pressure
    Qs = Cd*crack_area*(2*abs(stack_pressure_diff)/rho_i)**0.5 # Flow rate due to stack pressure
    Qt = (Qw**2 + Qs**2)**0.5 # Total flow rate
 
    total_air_flow = Qt*crack_length_total/crack_length 
    R_a_min = total_air_flow/V 

    # Ventilation
    DeltaT_vent = T_i - T_sp_vent
    comp_dtv_low = DeltaT_vent > 0 and DeltaT_vent < 4
    comp_dtv_high = DeltaT_vent >= 4
    R_a = R_a_min + comp_dtv_low*(R_a_max - R_a_min)/4*DeltaT_vent + comp_dtv_high*(R_a_max-R_a_min)

    QV_i_e = R_a*V*rho_i*c_i*(T_i - T_ext) # Internal air to outside air [J/s]
    QP_i_e = R_a*V*H_fg*(C_w - Cw_ext) # Latent heat loss due to leakiness
 
    MW_i_e = R_a*(C_w - Cw_ext)

    ##      Solar radiation
    # We first define the solar elevation angle that determines that absorption of solar radiation. Notation: r is direct radiation, f is diffuse radiation, whilst VIS and NIR stand for visible and near infra-red respectively.

    gamma = np.deg2rad(360.*(day(t) -  80.)/365.) # Year angle [rad] --- day counts from January 1st
    eqn_time = -7.13*np.cos(gamma) - 1.84*np.sin(gamma) - 0.69*np.cos(2.* gamma) + 9.92*np.sin(2.*gamma) # Equation of time [min]
    az = np.deg2rad(360.*((t/(3600.)%24.) + eqn_time/60. - 12.)/24.) # Azimuth [rad]
    delta = np.deg2rad(0.38 - 0.77*np.cos(gamma) + 23.27*np.cos(gamma)) # Declination angle [rad]
    lat = np.deg2rad(latitude)
    angler = np.arcsin(np.sin(lat)*np.sin(delta) + np.cos(lat)*np.cos(delta)*np.cos(az)) # Angle of elevation [rad]
    angle = np.rad2deg(angler)

    # Radiation from artifical lighting
    QS_al_NIR = 0. # no artificial lighting
    QS_al_VIS = 0. 

    # Solar radiation incident on the cover
    QS_tot_rNIR = 0.5*SurfaceArea@climate[n, 4:12] # Direct 
    QS_tot_rVIS = 0.5*SurfaceArea@climate[n, 4:12]
    QS_tot_fNIR = 0.5*SurfaceArea@climate[n, 12:20] # Diffuse
    QS_tot_fVIS = 0.5*SurfaceArea@climate[n, 12:20]

    # Transmitted solar radiation
    QS_int_rNIR = tau_c_NIR*QS_tot_rNIR # J/s total inside greenhouse
    QS_int_rVIS = tau_c_VIS*QS_tot_rVIS
    QS_int_fNIR = tau_c_NIR*QS_tot_fNIR
    QS_int_fVIS = tau_c_VIS*QS_tot_fVIS 

    # Solar radiation absorbed by the cover and the obstructions
    QS_c = alph_c*(QS_tot_rNIR + QS_tot_rVIS + QS_tot_fNIR + QS_tot_fVIS) # J/s
    QS_i = a_obs*(QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS) 

    # Solar radiation absorbed by the vegetation
    # Area = A_v i.e. planted area
    # factor QS by A_v/A_f

    k_fNIR = 0.27 # Near-IR diffuse extinction coefficient [-]
    a_v_fNIR = 0.65 - 0.65*np.exp(-k_fNIR*LAI) # Near-IR diffuse absorption coefficient [-]

    k_fVIS = 0.85 # Visible diffuse extinction coefficient [-]
    a_v_fVIS = 0.95 - 0.9*np.exp(-k_fVIS*LAI) # Visible diffuse absorption coefficient [-]

    k_rNIR = 0.25 + 0.38*np.exp(-0.12*angle) # Near-IR direct extinction coefficient [-]
    a_v_rNIR = 0.67 - 0.06*np.exp(-0.08*angle) - (0.68 - 0.5*np.exp(-0.11*angle))*np.exp(-k_rNIR*LAI) # Near-IR direct absorption coefficient [-]

    k_rVIS = 0.88 + 2.6*np.exp(-0.18*angle) # Visible direct extinction coefficient [-]
    a_v_rVIS = 0.94 - 0.95*np.exp(-k_rVIS*LAI) # Visible direct absorption coefficient [-]

    QS_v_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_v_rNIR*A_v/A_f
    QS_v_fNIR = (QS_int_fNIR*(1 - a_obs))*a_v_fNIR*A_v/A_f
    QS_v_NIR = (QS_v_rNIR + QS_v_fNIR) # factor as planted area not entire floor

    QS_v_rVIS = (QS_int_rVIS*(1 - a_obs) + QS_al_VIS)*a_v_rVIS*A_v/A_f
    QS_v_fVIS = (QS_int_fVIS*(1 - a_obs))*a_v_fVIS*A_v/A_f
    QS_v_VIS = (QS_v_rVIS + QS_v_fVIS) # Used for photosynthesis calc

    # Solar radiation absorbed by the mat
    a_m_fNIR = 0.05 + 0.91*np.exp(-0.5*LAI) # Near-IR diffuse absorption coefficient [-]
    a_m_fVIS = np.exp(-0.92*LAI) # Visible diffuse absorption coefficient [-]
    a_m_rNIR = 0.05 + 0.06*np.exp(-0.08*angle) + (0.92 - 0.53*np.exp(-0.18*angle))*np.exp(-(0.48 + 0.54*np.exp(-0.13*angle))*LAI) # Near-IR direct absorption coefficient [-]
    a_m_rVIS = np.exp(-(0.9 + 0.83*np.exp(-0.12*angle))*LAI) # Visible direct absorption coefficient [-]

    QS_m_rNIR = (QS_int_rNIR*(1 - a_obs) + QS_al_NIR)*a_m_rNIR*A_v/A_f
    QS_m_fNIR = QS_int_fNIR*(1 - a_obs)*a_m_fNIR*A_v/A_f # W
    QS_m_NIR = (QS_m_rNIR + QS_m_fNIR)

    QS_m_rVIS = (QS_int_rVIS*(1 - a_obs) + QS_al_VIS)*a_m_rVIS*A_v/A_f 
    QS_m_fVIS = QS_int_fVIS*(1 - a_obs)*a_m_fVIS*A_v/A_f
    QS_m_VIS = (QS_m_rVIS + QS_m_fVIS)

    # Solar radiation absorbed by the floor
    # factor by (A_f-A_v)/A_f

    QS_s_rNIR = QS_int_rNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_fNIR = QS_int_fNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_NIR = QS_s_rNIR + QS_s_fNIR

    QS_s_rVIS = QS_int_rVIS*(1-a_obs)*alphS_s*(A_f-A_v)/A_f 
    QS_s_fVIS = QS_int_fVIS*(1-a_obs)*alphS_s*(A_f-A_v)/A_f
    QS_s_VIS = QS_s_rVIS + QS_s_fVIS

    ## Transpiration
    QS_int = (QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS)*(1-a_obs)*A_v/A_f # J/s

    #  Vapour pressure deficit at leaf surface
    xa = C_w/rho_i #[-]
    xv = sat_conc(T_v)/rho_i #[-]
    vpd = atm*(xv/(xv + 0.622) - xa/(xa + 0.622)) # [Pa]

    # Stomatal resistance according to Stanghellini
    x = np.exp(-0.24*LAI) # [-]
    a_v_short = 0.83*(1 - 0.70*x)*(1 + 0.58*x**2)*(0.88 - x**2 + 0.12*x**(8/3)) # [-]Absorption for shortwave radiation
    I_s_bar = QS_int*a_v_short/(2*LAI) # [J/s] Mean radiation interacting with leaf surface

    Heavy_CO2 = I_s_bar > 0.
    r_i_CO2 = 1 + Heavy_CO2*6.1e-7*(C_c_ppm - 200)**2
    Heavy_vpd = vpd/1000 < 0.8
    r_i_vpd = Heavy_vpd*(1 + 4.3*(vpd/1000)**2) + (1 - Heavy_vpd)*3.8
    r_st = 82*((QS_int + 4.3)/(QS_int + 0.54))*(1 + 0.023*(T_v - T_k - 24.5)**2)*r_i_CO2*r_i_vpd #[s/m]

    hL_v_i = 2*LAI*H_fg/(rho_i*c_i)*(Le**(2/3)/HV + r_st/(rho_i*c_i))**(-1)

    QT_St = A_v*hL_v_i*(sat_conc(T_v) - C_w) # J/s

    QT_v_i = max(QT_St,0)
    
    ## Dehumidification
    MW_cc_i = 0 # No dehumidification included

    # CO2 exchange with outside
    MC_i_e = (R_a*(C_c - C_ce)) # [kg/m^3/s]

    day_hour_c=(hour/24-np.floor(hour/24))*24
    track=day_hour_c>6 and day_hour_c<20
    Value=added_CO2/Nz/3600./V

    MC_cc_i=Value*track

    ## Photosynthesis model - Vanthoor

    # Consider photosynthetically active radiation to be visible radiation

    T_25 = T_k + 25. # K

    I_VIS=QS_v_VIS # J/s incident on planted area

    PAR = I_VIS/heat_phot/N_A/A_v 

    # The number of moles of photosynthetically active photons per unit area of planted floor [mol{phot}/m^2/s]
    #J/s/(J/photon)/(photons/mol)/m^2 cf Vanthoor 2.3mumol(photons)/J

    Gamma = max((c_Gamma*(T_v - T_k)/LAI + 20*c_Gamma*(1 - 1/LAI)),0) # The CO2 compensation point [mol{CO2}/mol{air}]
    k_switch = C_buf_max # kg/m^2/s
    h_airbuf_buf = 1/(1 + np.exp(s_airbuf_buf*(C_buf - k_switch)))

    C_c_molar=(C_c/rho_i)*(M_a/M_c)
    C_stom = eta*C_c_molar # Stomatal CO2 concentration [mol{CO2}/mol{air}] 

    J_pot = LAI*J_max_25*np.exp(E_j*(T_v - T_25)/(R*T_v*T_25))*(1 + np.exp((S*T_25 - HH)/(R*T_25)))/(1 + np.exp((S*T_v - HH)/(R*T_v))) # [mol{e}/m^2{floor}s]
    J = (J_pot + alph*PAR - ((J_pot + alph*PAR)**2 - 4*theta*J_pot*alph*PAR)**0.5)/(2*theta)
    P = J*(C_stom - Gamma)/(4*(C_stom + 2*Gamma)) # Photosynthesis rate [mol{CO2}/s]
    Resp = P*Gamma/C_stom # Photorespiration rate

    MC_i_buf = (M_carb*h_airbuf_buf*(P - Resp)) # The net photosynthesis rate [kg{CH2O}/m^2/s]

    ## Crop growth model

    # Flow of carbohydrates from buffer to fruit, leaves and stem
    C_buf_min = 0.05*C_buf_max
    h_buforg_buf =1/(1 + np.exp(s_buforg_buf*(C_buf - C_buf_min)))

    # inhibition terms need temperatures in oC
    h_T_v = 1/(1 + np.exp(s_min_T*((T_v-T_k) - T_min_v)))/(1 + np.exp(s_max_T*((T_v-T_k) - T_max_v))) 
    h_T_v24 = 1/(1 + np.exp(s_min_T24*((T_vmean-T_k) - T_min_v24)))/(1 + np.exp(s_max_T24*((T_vmean-T_k) - T_max_v24))) 

    h_T_vsum = 0.5*(T_vsum/T_sum_end + ((T_vsum/T_sum_end)**2 + 1e-4)**0.5) - 0.5*(((T_vsum - T_sum_end)/T_sum_end)+(((T_vsum - T_sum_end)/T_sum_end)**2 + 1e-4)**0.5) 

    g_T_v24 = 0.047*(T_vmean - T_k) + 0.06

    MC_buf_fruit = (h_buforg_buf*h_T_v*h_T_v24*h_T_vsum*g_T_v24*rg_fruit) 
    MC_buf_leaf = (h_buforg_buf*h_T_v24*g_T_v24*rg_leaf)
    MC_buf_stem = (h_buforg_buf*h_T_v24*g_T_v24*rg_stem)

    # Growth respiration, which is CO2 leaving the buffer
    MC_buf_i = c_fruit_g*MC_buf_fruit + c_leaf_g*MC_buf_leaf + c_stem_g*MC_buf_stem
    #MC_buf_i = 0

    # Maintenance respiration
    MC_fruit_i = (c_fruit_m*Q_10**(0.1*(T_vmean - T_25))*C_fruit*(1 - np.exp(-c_RGR*R_fruit)))
    MC_leaf_i = (c_leaf_m*Q_10**(0.1*(T_vmean - T_25))*C_leaf*(1 - np.exp(-c_RGR*R_leaf)))
    MC_stem_i = (c_stem_m*Q_10**(0.1*(T_vmean - T_25))*C_stem*(1 - np.exp(-c_RGR*R_stem)))

    C_max_leaf = LAI_max/SLA
    MC_leaf_prune = max(C_leaf - C_max_leaf, 0)

    ## ODE equations
    
    # Temperature components
    dT_c_dt = (1/(A_c*cd_c))*(QV_i_c + QP_i_c - QR_c_f - QR_c_v - QR_c_m + QV_e_c - QR_c_sk + QS_c)

    dT_i_dt = (1/(V*rho_i*c_i))*(-QV_i_m - QV_i_v - QV_i_f - QV_i_c - QV_i_e - QV_i_p + QS_i)

    dT_v_dt = (1/(c_v*A_v*msd_v))*(QV_i_v - QR_v_c - QR_v_m - QR_v_p + QS_v_NIR - QT_v_i)

    dT_m_dt = (1/(A_m*c_m))*(QV_i_m + QP_i_m - QR_m_v - QR_m_c - QR_m_p - QD_m_p + QS_m_NIR)
    dT_p_dt = (1/(A_p*c_p))*(QD_m_p + QV_i_p + QP_i_p - QR_p_f - QR_p_v - QR_p_m)
    dT_f_dt = (1/(rhod_s[0]*A_f*c_s[0]*l_s[0]))*(QV_i_f + QP_i_f - QR_f_c - QR_f_p - QD_sf1 + QS_s_NIR)
    dT_s1_dt = (1/(rhod_s[1]*c_s[1]*l_s[1]*A_f))*(QD_sf1-QD_s12)
    dT_s2_dt = (1/(rhod_s[2]*c_s[2]*l_s[2]*A_f))*(QD_s12-QD_s23)
    dT_s3_dt = (1/(rhod_s[3]*c_s[3]*l_s[3]*A_f))*(QD_s23-QD_s34)
    dT_s4_dt = (1/(rhod_s[4]*c_s[4]*l_s[4]*A_f))*(QD_s34-QD_s45)

    # Water vapour
    dC_w_dt = (1/(V*H_fg))*(QT_v_i-QP_i_c-QP_i_f-QP_i_m-QP_i_p) - MW_i_e + MW_cc_i
    #dC_wdt = -MW_i_e

    # Carbon Dioxide
    dC_c_dt = MC_cc_i - MC_i_e + (M_c/M_carb)*(A_v/V)*(MC_buf_i + MC_fruit_i + MC_leaf_i + MC_stem_i - MC_i_buf)

    # Plant growth control
    dT_vmean_dt = 1/86400*(T_v - T_vmean) 
    dT_vsum_dt = 1/86400*(T_v - T_k) 

    # Plant carbon exchange
    dC_buf_dt = MC_i_buf - MC_buf_fruit - MC_buf_leaf - MC_buf_stem - MC_buf_i
    dC_fruit_dt = MC_buf_fruit - MC_fruit_i
    dC_leaf_dt = MC_buf_leaf - MC_leaf_i - MC_leaf_prune
    dC_stem_dt = MC_buf_stem - MC_stem_i

    # Plant growth
    dR_fruit_dt = (dC_fruit_dt/C_fruit - R_fruit)
    dR_leaf_dt = ((dC_leaf_dt + MC_leaf_prune)/C_leaf - R_leaf)
    dR_stem_dt = (dC_stem_dt/C_stem - R_stem)
    
    return np.array([dT_c_dt,dT_i_dt,dT_v_dt,dT_m_dt,dT_p_dt,dT_f_dt,dT_s1_dt,
                     dT_s2_dt,dT_s3_dt,dT_s4_dt,dT_vmean_dt,dT_vsum_dt,dC_w_dt,
                     dC_c_dt,dC_buf_dt,dC_fruit_dt,dC_leaf_dt,dC_stem_dt,
                     dR_fruit_dt,dR_leaf_dt,dR_stem_dt])
