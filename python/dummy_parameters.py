# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 10:35:04 2021

@author: rmw61
"""

import numpy as np

# Constants
g = 9.81
nu = 15.1e-6
ias = 0.4
lam = 0.025
cd_c = 7932.5
R = 8.314
atm = 1.013e5
M_w = 0.018
M_a = 0.029
H_fg = 2437000.
c_i = 1003.2
T_k = 273.15
N_A = 6.02214e+23
heat_phot = 3.6368e-19

# Geometry
A_c = 150.
A_f = 100.
A_p = 100.
A_m = A_p
V = 200.
l_f = 0.1
AF_g = 0.5
A_v = AF_g*A_p
A_l = 50.
l_m = 0.01

# Vegetation
c_v = 3500.
msd_v = 1.5
LAI = 0.5

# Critical dimensions
d_c = 2.0
d_f = 1.3
d_v = 0.1
d_m = 0.1
d_p = 1.

# View Factors
F_c_f = 0.6
F_f_c = 0.6
F_c_v = 0.2
F_c_m = 0.2
F_c_l = 0.
F_l_c = 0.1
F_l_v = 0.3
F_l_m = 0.3
F_l_p = 0.3
F_p_f = 0.3
F_p_v = 0.4
F_p_m = 0.4
F_p_l = 0.
F_v_c = 0.1
F_v_p = 0.5
F_v_m = 0.5
F_v_l = 0.
F_m_c = 0.25
F_m_p = 0.75
F_m_v = 0.5
F_m_l = 0.
F_f_p = 0.3

# Emissivity/Reflectivity
eps_c = 0.9
eps_f = 0.9
eps_v = 0.9
eps_m = 0.9
eps_p = 0.9
eps_l = 1.
rho_c = 0.1
rho_f = 0.1
rho_v = 0.1
rho_m = 0.1
rho_p = 0.1
rho_l = 0.

# Conduction
lam_c = [0.05, 50.0, 0.5, 1.5, 1.5, 1.5]
lam_f = 0.5
lam_p = 0.25
l_c = [0.005, 0.02, 0.05, 0.2, 0.3, 0.6]
c_c = [1000., 500., 800., 1200., 1200., 1200.]
c_f = 302400.
c_p = 12525.
rhod_c = [1000., 8000., 2400., 1800., 1800., 1800.]
#rho_fl = 2400

T_ss = 14. + T_k
T_fl = 18. + T_k

# Ventilation
ACH = 0.001

# Lights
P_al = 15000.
f_heat = 0.5
f_light = 1-f_heat
P_ambient_al = 500.
T_al = 25.+T_k

# Moisture
dsat = 0.5
c_m = 9000
dehumidify = 0.02


