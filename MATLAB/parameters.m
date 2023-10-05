%--------------------
%------ INPUTS ------
%--------------------

% CALCULATION CONTROL
% ===================
deltaT = 60; % time step size [s]
sim_days = 364; % number of days to simulate [-] NB 2016 Leap Year
Nz = 1; % number of zones [-]

% CONSTANTS
% =========
sigm = 5.67e-8; % Stefan-Boltzmann constant [W/m^2/K^4]
T_k = 273.15; % zero celsius [K]
g = 9.81; % acceleration due to gravity [m/s^2]
atm = 1.013e5; % standard atmospheric pressure [Pa]
latitude = 53.193583; % latitude of greenhouse
longitude = -2.434920; % longitude of greenhouse
N_A = 6.02214e+23; % Avogadro's number
M_a = 0.029; % molar mass of dry air [kg/mol]                                                                 
lambda = 0.025; % thermal conductivity of air [W/m/K]
c_i = 1003.2; % heat capacity of humid air [J/kg/K]                        
H_fg = 2437000; % latent heat of condensation of water [J/kg]
Le = 0.819; % Lewis number [-] 
R = 8.314; % gas constant [J/mol/K]
M_w = 0.018; % molar mass of water [kg/mol]
M_c = 0.044; % molar mass of CO2 [kg/mol]
M_carb = 0.03; % molar mass of CH2O [kg/mol]
nu = 15.1e-6; % kinematic viscosity [m^2/s]
rho_w = 1000; % density of water [kg/m^3]

% GEOMETRY
% ========
A_f = 250; % greenhouse floor area [m^2]
V = 1000; % greenhouse volume [m^3]
% surface areas NE Wall NE Roof SE Wall SE Roof SW Wall SW Roof NW Wall NW Roof [m^2]
SurfaceArea = [40 0 75 135 40 0 75 135]; 
A_c = sum(SurfaceArea, 2); % cover areas [m^2]
a_obs = 0.05; % fraction of solar radiation hitting obstructions [-]   
H = 5; % Height of greenhouse

% AIR CHARACTERISTICS
% ===================
int_air_speed = 0.5; % internal air speed [m/s]              
R_a_max = 30/3600; % ventilation air change rate [1/s] 
T_sp_vent = 25 + T_k; % Ventilation set-point

% COVER
% =====================
% Glass
eps_ce = 0.85; % far-IR emissivity, outer surface [-]
eps_ci = 0.85; % far-IR emissivity, inner surface [-]
tau_c = 0.0; % far-IR transmissivity (0.0) [-]
rho_ci = 0.15; % far-IR reflectivity, inner surface (0.1) [-]  
alph_c = 0.04; % solar absorptivity, taking 'perpendicular' values [-]
tau_c_NIR = 0.85; % near-IR transmissivity of cover (0.84) [-]  
tau_c_VIS = 0.85; % visible transmissivity of cover [-] 
d_c = 1.5; % characteristic length of cover [m]
cd_c = 8736; % cover heat capacity per unit area [J/m^2/K]

% FLOOR 
% =====================
lam_s=[1.7 0.85 0.85 0.85 0.85]; %thermal conductivity of soil layers [W/mK] Concrete, Soil, Clay
c_s=[880 1081 1081 1081 1081]; %specific heat of soil layers [J/kgK]
l_s=[0.02 0.05 0.1 0.25 1.0]; %thickness of soil layers [m]
rhod_s=[2300 1500 1600 1600 1600]; %density of soil layers [kg/m^3]
rho_s=0.85; %far-IR reflectance of floor [-]
eps_s=0.95; %far-IR emmittance of floor [-]
rhoS_s=0.5; %solar reflectance of floor [-]
alphS_s=0.5; %solar absorptance of floor [-]
d_f=0.5; %characteristic floor length [m]
T_ss = 14.0 + T_k; %deep soil temperature [K]

% VEGETATION CHARACTERISTICS
% ==========================
c_v = 4180; % heat capacity of vegetation [J/kgK]
k_l = 0.94; % long-wave extinction coefficient [-]
rho_v = 0.22; % far-IR reflectivity of vegetation [-]
eps_v = 0.95; % far-IR emissivity of vegetation [-]
rhoS_v=0.35; %solar reflectance of vegetation [-]
d_v = 0.1; % characteristic leaf length [m]
p_v = 0.75; % cultivated fraction of floor 
msd_v = 1.326; % surface density [kg/m^2]

% TRAY/MAT CHARACTERISTICS
% ========================
A_p = p_v*A_f; % Area of cultivated floor [m^2]
A_v = A_p; % Area of plants [m^2]
A_m = A_p; % Area of mat for conduction to tray [m^2]
d_p = 1; % characteristic dimension of tray (width)
d_m = 0.1; % characteristic dimension of mat (width)
lam_m = 0.5; % thermal conductivity of mat [W/mK]
lam_p = 0.2; % thermal conductivity of plastic tray [W/mK]
c_m = 45050; % specific heat of mat assumed 25% saturated [J/m^2K]  
c_p = 10020; % specific heat of tray [J/m^2K]
l_m = 0.03; % thickness of mat [m]
l_p = 0.005; % thickness of tray [m]
rhod_m = 720; % density of mat [kg/m^3]
rhod_p = 1200; % density of tray [kg/m^3]
rho_m = 0.05; % far-IR reflectivity of mat [-]
rho_p = 0.05; % far-IR reflectivity of tray
eps_m = 0.95; % far-IR emissivity of mat [-]
eps_p = 0.95; % far-IR emissivity of tray

% PHOTOSYNTHESIS MODEL - VANTHOOR
% ===============================
c_Gamma = 1.7e-6; % effect of canopy temperature on CO2 compensation point [mol{CO2}/mol{air}/K]
J_max_25 = 210e-6; % maximum rate of electron transport at 25 C [mol{e}/m^2{leaf}/s]
alph = 0.385; % conversion factor from photons to electrons [mol{e}/mol{phot}]
C_buf_max = 0.02; % maximum buffer capacity per unit area of cultivated floor [kg{CH2O}/m^2/s]
theta = 0.7; % degree of curvature of the electron transport rate [-]
S = 710; % entropy term for J_pot calculation [J/mol/K]
HH = 22e4; % deactivation energy for J_pot calculation [J/mol]
E_j = 37e3; % activation energy for J_pot calculation [J/mol]
heat_phot = 3.6368e-19; % conversion rate from incident energy to number of photons [num{photons}/J]
eta = 0.67; % conversion factor from CO2 in the air to CO2 in the stomata [-]
s_airbuf_buf = 5e2; % differential switch function slope for maximum buffer capacity [m^2/kg]
s_buforg_buf = -5e3; % differential switch function slope for minimum buffer capacity [m^2/kg]
s_min_T = -0.8690; % differential switch function slope for minimum photosynthesis instantaneous temperature [1/degC]
s_max_T = 0.5793; % differential switch function slope for maximum photosynthesis instantaneous temperature [1/degC]
s_min_T24 = -1.1587; % differential switch function slope for minimum photosynthesis mean 24 hour temperature [1/degC]
s_max_T24 = 1.3904; % differential switch function slope for maximum photosynthesis mean 24 hour temperature [1/degC]
s_prune = -50; % differential switch function slope for leaf pruning [m^2/kg]

% CROP GROWTH MODEL
% ===================
added_CO2 = 0; % mass of CO2 pumped in per hour [kg/h] (100)
SLA = 26.6; % specific leaf area index [m^2{leaf}/kg{CH2O}]
LAI_max = 5; % the maximum allowed leaf area index [m^2{leaf}/m^2{floor}]
Q_10 = 2; % see parameters for de Zwart model above [-]
rg_fruit = 0.328e-6; % potential fruit growth rate coefficient at 20 deg C [kg{CH2O}/m^2/s]
rg_leaf = 0.095e-6; % potential leaf growth rate coefficient at 20 deg C [kg{CH2O}/m^2/s]
rg_stem = 0.074e-6; % potential stem growth rate coefficient at 20 deg C [kg{CH2O}/m^2/s]
c_fruit_g = 0.27; % fruit growth respiration coefficient [-]
c_fruit_m = 1.16e-7; % fruit maintenance respiration coefficient [1/s]
c_leaf_g = 0.28; % leaf growth respiration coefficient [-]
c_leaf_m = 3.47e-7; % leaf maintenance respiration coefficient [1/s]
c_stem_g = 0.30; % stem growth respiration coefficient [-]
c_stem_m = 1.47e-7; % stem maintenance respiration coefficient [1/s]
c_RGR = 2.85e6; % regression coefficient in maintenance respiration function [s]
T_min_v24 = 12 ; % between base temperature and first optimal temperature for 24 hour mean [oC]
T_max_v24 = 27 ; % between second optimal temperature and maximum temperature for 24 hour mean [oC]
T_min_v = 6 ; % between base temperature and first optimal temperature [oC]
T_max_v = 40 ; % between second optimal temperature and maximum temperature [oC]
T_sum_end = 1035; % the temperature sum at which point the fruit growth rate is maximal [oC]


 % PARAMETERS FOR INFILTRATION CALCULATION
 % =======================================
 c = 0.35; % terrain factor, see Awbi, Chapter 3, Table 3.2
 a = 0.25; % terrain factor, see Awbi, Chapter 3, Table 3.2
 Cp = 0.62; % static pressure coefficient - for wind perpendicular to gap
 Cd = 0.61; % sharp edge orifice, see Awbi
 crack_length = 1; % typical estimate
 crack_width = 0.001; % typical estimate
 crack_area = crack_length*crack_width;
 crack_length_total = 350; % 