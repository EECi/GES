 function dT = derivatives(t, T)
%% IMPORT DATA, OUTPUT DAY NUMBER

global climate;
persistent day;

%% Read in and initialise parameters for model

parameters;

% Time control parameters
n = ceil(t/deltaT);
if t == 0
    n = 1;
    day = 1;
 end

day_new = ceil(t/(86400)); % Day number
if day_new == day + 1
    disp(strcat(num2str(day_new), '/', num2str(sim_days)));
end
day = day_new;

hour = floor(t/3600) + 1;



% Initial Values
T_c = T(1, 1); % Cover temperature [K]
T_i = T(2, 1); % Internal temperature [K]
T_v = T(3, 1); % Vegetation temperature [K]
T_m = T(4, 1); % Temperature of mat [K]
T_p = T(5, 1); %Temperature of tray [K] 
T_f = T(6, 1); % Temperature of floor [K]
T_s1 = T(7, 1); % Temperature of first layer of soil [K]
T_s2 = T(8, 1); % Temperature of second layer of soil [K]
T_s3 = T(9, 1); % Temperature of third layer of soil [K]
T_s4 = T(10, 1); % Temperature of fourth layer of soil [K]
%
T_vmean = T(11, 1); % 24 hour mean vegetation temperature [K]
T_vsum = T(12, 1); % Vegetation temperature sum [degC]                        
%
C_w = T(13, 1); % Density of water vapour [kg/m^3]
C_c = T(14, 1); % Density of carbon dioxide [kg/m^3]
C_buf = T(15, 1); % Mass of carbohydrate in buffer per unit per unit area of cultivated floor [kg/m^2]
C_fruit = T(16, 1); % Mass of carbohydrate in fruit per unit per unit area of cultivated floor [kg/m^2]
C_leaf = T(17, 1); % Mass of carbohydrate in leaves per unit per unit area of cultivated floor [kg/m^2]
C_stem = T(18, 1); % Mass of carbohydrate in stem per unit per unit area of cultivated floor [kg/m^2]
R_fruit = T(19, 1); % Relative growth rate of fruit averaged over 5 days [1/s]
R_leaf = T(20, 1); % Relative growth rate of leaf averaged over 5 days [1/s]
R_stem = T(21, 1); % Relative growth rate of stem averaged over 5 days [1/s]
%
T_e = climate(n, 1) + T_k; % External temperature [K]
T_sk = climate(n, 2) + T_k; % Sky temperature [K]
wind_speed = climate(n, 3)*ones(Nz, 1); % External wind speed [m/s]
rel_humid = climate(n, 4)/100; % External relative humidity [-]
C_we = rel_humid.*sat_conc(T_e); % External density of water vapour [kg/m^3]
p_w = C_w*R.*T_i/M_w; % Partial pressure of water [Pa]
rho_i = ((atm - p_w)*M_a + p_w*M_w)./(R*T_i); % Density of air [kg/m^3]   
LAI = SLA.*C_leaf; % Leaf area index
C_ce = 4.0e-4*M_c*atm/(R*T_e); % External carbon dioxide concentration [kg/m^3]
C_c_ppm = C_c*R.*T_i/(M_c*atm)*1e6; % CO2 concentration in ppm

%% HEAT EXCHANGE TERMS

% Convection including evaporation/condensation
% Radiation
% Conduction
% Transpiration

%%      Convection

% The Grashof and Reynolds numbers are calculated, and then the Nusselt
% number calculated as the larger of the 'free' Nusselt number, based on
% the Grashof number, and the 'forced' Nusselt number, based on the
% Reynolds number. The heat transferred per unit time is then calculated
% based on this. The first subscript is where heat is transferred from, and
% the second subscript where heat is transferred to.

% Cover to external air
Gr = g*d_c^3/(T_e*nu^2)*abs(T_c - T_e);
Re = wind_speed*d_c/nu;
Nu = lamorturb(Gr, Re);
QV_c_e = A_c.*Nu*lambda.*(T_c - T_e)/d_c;

% Internal air to cover
Gr = g*d_c^3./(T_i*nu^2).*abs(T_i - T_c);
Re = int_air_speed*d_c/nu;
[Nu, Sh] = lamorturb(Gr, Re, Le);
QV_i_c = A_c.*Nu*lambda.*(T_i - T_c)/d_c;
QP_i_c = A_c*H_fg./(rho_i*c_i).*Sh/Le*lambda/d_c.*(C_w - sat_conc(T_c));
QP_i_c = max(QP_i_c,0); % assumed no evaporation from the cover, only condensation

%Internal air to floor
Gr = g*d_f^3./(T_i*nu^2).*abs(T_i - T_f); 
Re = int_air_speed*d_f/nu;
[Nu, Sh] = lamorturb(Gr, Re, Le);
QV_i_f = A_f.*Nu*lambda.*(T_i - T_f)/d_f;
QP_i_f = A_f*H_fg./(rho_i*c_i).*Sh/Le*lambda/d_f.*(C_w - sat_conc(T_f)); 
QP_i_f = max(QP_i_f,0); % assumed no evaporation from the floor, only condensation

% Mat to internal air
Gr = g*d_m^3./(T_i*nu^2).*abs(T_m - T_i);
Re = int_air_speed*d_m/nu;
[Nu, Sh] = lamorturb(Gr, Re, Le);
A_m_wool = 0.75*A_m; % Area of mat exposed
A_m_water = 0.25*A_m; % assumed 25% saturated
QV_m_i = A_m_wool.*Nu*lambda.*(T_m - T_i)/d_m;
QP_m_i = A_m_water*H_fg./(rho_i*c_i).*Sh/Le*lambda/d_m.*(sat_conc(T_m) - C_w);

% Tray to internal air
Gr = g*d_p^3./(T_i*nu^2).*abs(T_p - T_i);
Re = int_air_speed*d_p/nu;
[Nu, ~] = lamorturb(Gr, Re, Le);
QV_p_i = A_p.*Nu*lambda.*(T_p - T_i)/d_p;
QP_p_i = 0; % assumed no condensation/evaporation from tray

% Vegetation to internal air
Gr = g*d_v^3./(T_i*nu^2).*abs(T_v - T_i);
Re = int_air_speed*d_v/nu;
Nu = lamorturb(Gr, Re);
A_v_exp = LAI*A_v; % Area vegetation exposed
QV_v_i = A_v_exp.*Nu*lambda.*(T_v - T_i)/d_v;
HV = Nu*lambda/d_v; % Used in the transpiration term

%%      Far-IR radiation

% In the following equations the emissivity is used at times when the
% relevant quantity is the absorptivity, since these are equal according to
% Kirchoff's law. The constant k is the pre-factor that takes into account
% emission, reflection, transmission and absorption. The F's are view
% factors --- in the case of the cover they are calculated using
% reciprocity since the precise geometry of the cover is not known. Once
% again, the first subscript is the element that is radiating, whilst the
% second subscript is the element that is absorbing. 

% View factors
A_vvf=min(LAI*p_v*A_f,p_v*A_f);
A_c_roof = 271;

F_f_c = 1-p_v;  % Floor to cover
F_f_p = p_v; % Floor to tray
F_c_f = A_f/A_c_roof*F_f_c; % Cover to floor
F_c_v = min((1-F_c_f)*LAI,(1-F_c_f)); % Cover to vegetation
F_c_m = max((1-F_c_f)*(1-LAI),0); % Cover to mat
F_v_c = 0.5; % Vegetation to cover
F_v_m = 0.5; % Vegetation to mat
F_v_p = 0; % Vegetation to tray
F_m_c = max((1-LAI),0.0); % Mat to cover
F_m_v = 1-F_m_c; % Mat to vegetation
F_m_p = 0; % Mat to tray
F_p_v = 0; % Tray to vegetation
F_p_m = 0; % Tray to mat
F_p_f = 1.0; % Tray to floor

% Cover to sky
k = eps_ce;
QR_c_sk = k*sigm*A_c.*(T_c.^4 - T_sk^4); % J/s

% Cover to vegetation
k = eps_ci*eps_v./(1 - rho_ci*rho_v*F_c_v.*F_v_c);
QR_c_v = k*sigm.*A_c_roof.*F_c_v.*(T_c.^4 - T_v.^4);

% Cover to mat
k = eps_ci*eps_m./(1 - rho_ci*rho_m*F_c_m.*F_m_c);
QR_c_m = k*sigm.*A_c_roof.*F_c_m.*(T_c.^4 - T_m.^4);

% Cover to floor
k = eps_ci*eps_s./(1 - rho_ci*rho_s*F_c_f.*F_f_c);
QR_c_f = k*sigm.*A_c_roof.*F_c_f.*(T_c.^4 - T_f.^4);

% Vegetation to cover
k = eps_v*eps_ci./(1 - rho_v*rho_ci*F_v_c.*F_c_v);
QR_v_c = k*sigm.*A_vvf*2.*F_v_c.*(T_v.^4 - T_c.^4);

% Vegetation to mat
k = eps_v*eps_m./(1 - rho_v*rho_m*F_v_m.*F_m_v);
QR_v_m = k*sigm.*A_vvf*2.*F_v_m.*(T_v.^4 - T_m.^4);

% Vegetation to tray
k = eps_v*eps_p./(1 - rho_v*rho_p*F_v_p.*F_p_v);
QR_v_p = k*sigm.*A_vvf*2.*F_v_p.*(T_v.^4 - T_p.^4);

% Mat to cover
k = eps_m*eps_ci./(1 - rho_m*rho_ci*F_m_c.*F_c_m);
QR_m_c = k*sigm.*A_m.*F_m_c.*(T_m.^4 - T_c.^4);

% Mat to vegetation
k = eps_m*eps_v./(1 - rho_m*rho_v*F_m_v.*F_v_m);
QR_m_v = k*sigm.*A_m.*F_m_v.*(T_m.^4 - T_v.^4);

% Mat to tray
k = eps_m*eps_p./(1 - rho_m*rho_p*F_m_p.*F_p_m);
QR_m_p = k*sigm.*A_m.*F_m_p.*(T_m.^4 - T_p.^4);

% Tray to vegetation
k = eps_p*eps_v./(1 - rho_p*rho_v*F_p_v.*F_v_p);
QR_p_v = k*sigm.*A_p.*F_p_v.*(T_p.^4 - T_v.^4);

% Tray to mat
k = eps_p*eps_m./(1 - rho_p*rho_m*F_p_m.*F_m_p);
QR_p_m = k*sigm.*A_p.*F_p_m.*(T_p.^4 - T_m.^4);

% Tray to floor
k = eps_p*eps_s./(1 - rho_p*rho_s*F_p_f.*F_f_p);
QR_p_f = k*sigm.*A_p.*F_p_f.*(T_p.^4 - T_f.^4);

% Floor to cover
k = eps_s*eps_ci./(1 - rho_s*rho_ci*F_c_f.*F_f_c);
QR_f_c = k*sigm.*A_f.*F_f_c.*(T_f.^4 - T_c.^4); 

% Floor to tray
k = eps_s*eps_p./(1 - rho_s*rho_p*F_f_p.*F_p_f);
QR_f_p = k*sigm.*A_f.*F_f_p.*(T_f.^4 - T_p.^4); 


%%      Solar radiation

% We first define the solar elevation angle that determines that absorption
% of solar radiation. Notation: r is direct radiation, f is diffuse
% radiation, whilst VIS and NIR stand for visible and near infra-red
% respectively.

gamma = 360*(day -  80)/365; % Year angle [deg] --- day counts from January 1st
eqn_time = -7.13*cosd(gamma) - 1.84*sind(gamma) - ...
    0.69*cosd(2*gamma) + 9.92*sind(2*gamma); % Equation of time [min]
az = 360*(rem(t/(3600), 24) + eqn_time/60 - 12)/24; % Azimuth [deg]
delta = 0.38 - 0.77*cosd(gamma) + 23.27*cosd(gamma); % Declination angle [deg]
angle = asind(sind(latitude)*sind(delta) + ...
    cosd(latitude)*cosd(delta)*cosd(az)); % Angle of elevation [deg]

% Radiation from artifical lighting
QS_al_NIR = 0; % no artificial lighting
QS_al_VIS = 0;

% Solar radiation incident on the cover
QS_tot_rNIR = 0.5*SurfaceArea*climate(n, 5:12)'; % Direct (J/s)
QS_tot_rVIS = 0.5*SurfaceArea*climate(n, 5:12)';
QS_tot_fNIR = 0.5*SurfaceArea*climate(n, 13:20)'; % Diffuse
QS_tot_fVIS = 0.5*SurfaceArea*climate(n, 13:20)';

% Transmitted solar radiation
QS_int_rNIR = tau_c_NIR*QS_tot_rNIR; % J/s total inside greenhouse
QS_int_rVIS = tau_c_VIS*QS_tot_rVIS;
QS_int_fNIR = tau_c_NIR*QS_tot_fNIR;
QS_int_fVIS = tau_c_VIS*QS_tot_fVIS; % 

% Solar radiation absorbed by the cover and the obstructions
QS_c = alph_c*(QS_tot_rNIR + QS_tot_rVIS + QS_tot_fNIR + QS_tot_fVIS); % J/s
QS_i = a_obs*(QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS); 

% Solar radiation absorbed by the vegetation
% Area = A_v i.e. planted area
% factor QS by A_v/A_f

k_fNIR = 0.27; % Near-IR diffuse extinction coefficient [-]
a_v_fNIR = 0.65 - 0.65*exp(-k_fNIR*LAI); % Near-IR diffuse absorption coefficient [-]
%
k_fVIS = 0.85; % Visible diffuse extinction coefficient [-]
a_v_fVIS = 0.95 - 0.9*exp(-k_fVIS*LAI); % Visible diffuse absorption coefficient [-]
%
k_rNIR = 0.25 + 0.38*exp(-0.12*angle); % Near-IR direct extinction coefficient [-]
a_v_rNIR = 0.67 - 0.06*exp(-0.08*angle) - (0.68 - 0.5*exp(-0.11*angle))* ...
    exp(-k_rNIR*LAI); % Near-IR direct absorption coefficient [-]
%
k_rVIS = 0.88 + 2.6*exp(-0.18*angle); % Visible direct extinction coefficient [-]
a_v_rVIS = 0.94 - 0.95*exp(-k_rVIS*LAI); % Visible direct absorption coefficient [-]

QS_v_rNIR = (QS_int_rNIR.*(1 - a_obs) + QS_al_NIR)*a_v_rNIR*A_v/A_f;
QS_v_fNIR = (QS_int_fNIR.*(1 - a_obs))*a_v_fNIR*A_v/A_f;
QS_v_NIR = (QS_v_rNIR + QS_v_fNIR); % factor as planted area not entire floor

QS_v_rVIS = (QS_int_rVIS.*(1 - a_obs) + QS_al_VIS)*a_v_rVIS*A_v/A_f;
QS_v_fVIS = (QS_int_fVIS.*(1 - a_obs))*a_v_fVIS*A_v/A_f;
QS_v_VIS = (QS_v_rVIS + QS_v_fVIS); % Used for photosynthesis calc

% Solar radiation absorbed by the mat
a_m_fNIR = 0.05 + 0.91*exp(-0.5*LAI); % Near-IR diffuse absorption coefficient [-]
a_m_fVIS = exp(-0.92*LAI); % Visible diffuse absorption coefficient [-]
a_m_rNIR = 0.05 + 0.06*exp(-0.08*angle) + (0.92 - 0.53*exp(-0.18*angle))* ...
    exp(-(0.48 + 0.54*exp(-0.13*angle))*LAI); % Near-IR direct absorption coefficient [-]
a_m_rVIS = exp(-(0.9 + 0.83*exp(-0.12*angle))*LAI); % Visible direct absorption coefficient [-]

QS_m_rNIR = (QS_int_rNIR.*(1 - a_obs) + QS_al_NIR).*a_m_rNIR*A_v/A_f;
QS_m_fNIR = QS_int_fNIR.*(1 - a_obs).*a_m_fNIR*A_v/A_f; % W
QS_m_NIR = (QS_m_rNIR + QS_m_fNIR);

QS_m_rVIS = (QS_int_rVIS.*(1 - a_obs) + QS_al_VIS).*a_m_rVIS*A_v/A_f; 
QS_m_fVIS = QS_int_fVIS.*(1 - a_obs).*a_m_fVIS*A_v/A_f;
QS_m_VIS = (QS_m_rVIS + QS_m_fVIS);

% Solar radiation absorbed by the floor
% factor by (A_f-A_v)/A_f

QS_s_rNIR = QS_int_rNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f;
QS_s_fNIR = QS_int_fNIR*(1-a_obs)*alphS_s*(A_f-A_v)/A_f;
QS_s_NIR = QS_s_rNIR + QS_s_fNIR;

QS_s_rVIS = QS_int_rVIS*(1-a_obs)*alphS_s*(A_f-A_v)/A_f ;
QS_s_fVIS = QS_int_fVIS*(1-a_obs)*alphS_s*(A_f-A_v)/A_f;
QS_s_VIS = QS_s_rVIS + QS_s_fVIS;
%

%%      Transpiration

% Radiation
QS_int = (QS_int_rNIR + QS_int_rVIS + QS_int_fNIR + QS_int_fVIS)*(1-a_obs)*A_v/A_f; % J/s

% Vapour pressure deficit at leaf surface
xa = C_w./rho_i; %[-]
xv = sat_conc(T_v)./rho_i; %[-]
vpd = atm*(xv./(xv + 0.622) - xa./(xa + 0.622)); % [Pa]

% Stomatal resistance according to Stanghellini
x = exp(-0.24*LAI); % [-]
a_v_short = 0.83*(1 - 0.70*x).*(1 + 0.58*x.^2).*(0.88 - x.^2 + 0.12*x.^(8/3)); % [-]Absorption for shortwave radiation
I_s_bar = QS_int.*a_v_short./(2*LAI); % [J/s] Mean radiation interacting with leaf surface

Heavy_CO2 = I_s_bar > 0;
r_i_CO2 = 1 + Heavy_CO2*6.1e-7.*(C_c_ppm - 200).^2;
Heavy_vpd = vpd/1000 < 0.8;
r_i_vpd = Heavy_vpd.*(1 + 4.3*(vpd/1000).^2) + (1 - Heavy_vpd)*3.8;
r_st = 82*((QS_int + 4.3)./(QS_int + 0.54)).*(1 + 0.023*(T_v - T_k - 24.5).^2).*r_i_CO2.*r_i_vpd; %[s/m]

hL_v_i = 2*LAI*H_fg./(rho_i*c_i).*(Le^(2/3)./HV + r_st./(rho_i*c_i)).^(-1); %

QT_St = A_v.*hL_v_i.*(sat_conc(T_v) - C_w); % J/s

QT_v_i = max(QT_St,0);

%%      Conduction
% Floor layers 
% 
QD_sf1 = A_f*lam_s(1)/l_s(1).*(T_f - T_s1); 
QD_s12 = A_f*lam_s(2)/l_s(2).*(T_s1 - T_s2);
QD_s23 = A_f*lam_s(3)/l_s(3).*(T_s2 - T_s3);
QD_s34 = A_f*lam_s(4)/l_s(4).*(T_s3 - T_s4);
QD_s45 = A_f*lam_s(5)/l_s(5).*(T_s4 - T_ss);

% Mat to tray
QD_m_p = A_m*lam_p/l_m.*(T_m-T_p);

%%      Ventilation

% Leakage (equations for orifice flow from Awbi, Ventilation of Buildings, Chapter 3)
wind_speed_H = wind_speed*c*H^a; % Wind speed at height H
wind_pressure = Cp*0.5*rho_i.*wind_speed_H.^2; % Equals DeltaP for wind pressure
stack_pressure_diff = rho_i*g*H.*(T_i - T_e)./T_i; % DeltaP for stack pressure
% 
Qw = Cd*crack_area*(2*wind_pressure./rho_i).^0.5; % Flow rate due to wind pressure
Qs = Cd*crack_area*(2*abs(stack_pressure_diff)./rho_i).^0.5; % Flow rate due to stack pressure
Qt = (Qw.^2 + Qs.^2).^0.5; % Total flow rate
% 
total_air_flow = Qt.*crack_length_total/crack_length; % 
R_a_min = total_air_flow./V; % 

% Ventilation
DeltaT_vent = T_i - T_sp_vent;
comp_dtv_low = DeltaT_vent > 0 & DeltaT_vent < 4;
comp_dtv_high = DeltaT_vent >= 4;
R_a = R_a_min + comp_dtv_low.*(R_a_max - R_a_min)/4.*DeltaT_vent + comp_dtv_high.*(R_a_max-R_a_min);

QV_i_e = R_a.*V.*rho_i*c_i.*(T_i - T_e); % Internal air to outside air [J/s]
QP_i_e = R_a.*V*H_fg.*(C_w - C_we); % Latent heat loss due to leakiness


%% WATER VAPOUR EXCHANGE TERMS

MW_i_e = R_a.*(C_w - C_we);
MW_cc_i = 0; % Climate controller

%% CARBON DIOXIDE EXCHANGE TERMS

MC_i_e = (R_a.*(C_c - C_ce)); % [kg/m^3/s]

day_hour_c=(hour/24-floor(hour/24))*24;
track=day_hour_c>6 && day_hour_c<20;
Value=added_CO2/Nz/3600./V;

MC_cc_i=Value*track;


%% PHOTOSYNTHESIS MODEL - VANTHOOR

% Consider photosynthetically active radiation to be visible radiation

T_25 = T_k + 25; % K

I_VIS=QS_v_VIS; % J/s incident on planted area

PAR = I_VIS/heat_phot/N_A./A_v; 
% The number of moles of photosynthetically active photons per unit area of planted floor [mol{phot}/m^2/s]
%J/s/(J/photon)/(photons/mol)/m^2 cf Vanthoor 2.3mumol(photons)/J

Gamma = max((c_Gamma*(T_v - T_k)./LAI + 20*c_Gamma*(1 - 1./LAI)),0); % The CO2 compensation point [mol{CO2}/mol{air}]
k_switch = C_buf_max;% kg/m^2/s
h_airbuf_buf = 1./(1 + exp(s_airbuf_buf*(C_buf - k_switch)));

C_c_molar=(C_c./rho_i).*(M_a/M_c);
C_stom = eta*C_c_molar; % Stomatal CO2 concentration [mol{CO2}/mol{air}] 

J_pot = LAI*J_max_25.*exp(E_j*(T_v - T_25)./(R*T_v*T_25)).*...
    (1 + exp((S*T_25 - HH)/(R*T_25)))./(1 + exp((S*T_v - HH)./(R*T_v))); % [mol{e}/m^2{floor}s]
J = (J_pot + alph*PAR - sqrt((J_pot + alph*PAR).^2 - 4*theta*J_pot*alph.*PAR))/(2*theta);
P = J.*(C_stom - Gamma)./(4*(C_stom + 2*Gamma)); % Photosynthesis rate [mol{CO2}/s]
Resp = P.*Gamma./C_stom; % Photorespiration rate

MC_i_buf = (M_carb*h_airbuf_buf.*(P - Resp)); % The net photosynthesis rate [kg{CH2O}/m^2/s]


%% CROP GROWTH MODEL

% Flow of carbohydrates from buffer to fruit, leaves and stem
C_buf_min = 0.05*C_buf_max;
h_buforg_buf =1./(1 + exp(s_buforg_buf*(C_buf - C_buf_min)));

% inhibition terms need temperatures in oC
h_T_v = 1./(1 + exp(s_min_T*((T_v-T_k) - T_min_v)))./(1 + exp(s_max_T*((T_v-T_k) - T_max_v))); 
h_T_v24 = 1./(1 + exp(s_min_T24*((T_vmean-T_k) - T_min_v24)))./(1 + exp(s_max_T24*((T_vmean-T_k) - T_max_v24))); 

h_T_vsum = 0.5*(T_vsum./T_sum_end + sqrt((T_vsum./T_sum_end).^2 + 1e-4)) ...
    - 0.5*(((T_vsum - T_sum_end)/T_sum_end)+sqrt(((T_vsum - T_sum_end)/T_sum_end).^2 + 1e-4)); 

g_T_v24 = 0.047*(T_vmean - T_k) + 0.06;

MC_buf_fruit = (h_buforg_buf.*h_T_v.*h_T_v24.*h_T_vsum.*g_T_v24.*rg_fruit); 
MC_buf_leaf = (h_buforg_buf.*h_T_v24.*g_T_v24.*rg_leaf);
MC_buf_stem = (h_buforg_buf.*h_T_v24.*g_T_v24.*rg_stem);

% Growth respiration, which is CO2 leaving the buffer
MC_buf_i = c_fruit_g*MC_buf_fruit + c_leaf_g*MC_buf_leaf + c_stem_g*MC_buf_stem;

% Maintenance respiration
MC_fruit_i = (c_fruit_m*Q_10.^(0.1*(T_vmean - T_25)).*C_fruit.*(1 - exp(-c_RGR*R_fruit)));
MC_leaf_i = (c_leaf_m*Q_10.^(0.1*(T_vmean - T_25)).*C_leaf.*(1 - exp(-c_RGR*R_leaf)));
MC_stem_i = (c_stem_m*Q_10.^(0.1*(T_vmean - T_25)).*C_stem.*(1 - exp(-c_RGR*R_stem)));

C_max_leaf = LAI_max/SLA;
MC_leaf_prune = max(C_leaf - C_max_leaf, 0);

%% CALCULATING DERIVATIVES

% Temperature components
dT_c_dt = 1./(A_c*cd_c).*(QV_i_c + QP_i_c - QR_c_f - QR_c_v ...
    - QR_c_m - QV_c_e - QR_c_sk + QS_c); 
dT_i_dt = 1./(V.*rho_i*c_i).*(QV_m_i + QV_v_i - QV_i_f - QV_i_c - QV_i_e ...
    + QV_p_i + QS_i); 
dT_v_dt = 1./(c_v*A_v.*msd_v).*(-QV_v_i - QR_v_c - QR_v_m ...
    - QR_v_p - QT_v_i + QS_v_NIR); 
dT_m_dt = 1./(A_m*c_m).*(-QV_m_i - QP_m_i - QR_m_v - QR_m_c ...
    - QR_m_p - QD_m_p + QS_m_NIR);
dT_p_dt = 1./(A_p*c_p).*(QD_m_p - QV_p_i - QP_p_i - QR_p_f ...
    - QR_p_v - QR_p_m);
dT_f_dt = 1./(rhod_s(1)*A_f*c_s(1)*l_s(1)).*(QV_i_f + QP_i_f - ...
    QR_f_c - QR_f_p - QD_sf1 + QS_s_NIR);
dT_s1_dt = 1./(rhod_s(2)*c_s(2)*l_s(2)*A_f).*(QD_sf1 - QD_s12);
dT_s2_dt = 1./(rhod_s(3)*c_s(3)*l_s(3)*A_f).*(QD_s12 - QD_s23);
dT_s3_dt = 1./(rhod_s(4)*c_s(4)*l_s(4)*A_f).*(QD_s23 - QD_s34);
dT_s4_dt = 1./(rhod_s(5)*c_s(5)*l_s(5)*A_f).*(QD_s34 - QD_s45);

% Water Vapour
dC_w_dt = 1./(V*H_fg).*(QP_m_i + QP_p_i - QP_i_c - QP_i_f + QT_v_i) + ...
    MW_cc_i - MW_i_e;

% Carbon Dioxide
dC_c_dt = MC_cc_i - MC_i_e + (M_c/M_carb)*(A_v/V)*(MC_buf_i + MC_fruit_i + ...
    MC_leaf_i + MC_stem_i - MC_i_buf);

% Plant growth control
dT_vmean_dt = 1/86400*(T_v - T_vmean); 
dT_vsum_dt = 1/86400*(T_v - T_k); 

% Plant carbon exchange
dC_buf_dt = MC_i_buf - MC_buf_fruit - MC_buf_leaf - MC_buf_stem - MC_buf_i;
dC_fruit_dt = MC_buf_fruit - MC_fruit_i;
dC_leaf_dt = MC_buf_leaf - MC_leaf_i - MC_leaf_prune;
dC_stem_dt = MC_buf_stem - MC_stem_i;

% Plant growth
dR_fruit_dt = (dC_fruit_dt./C_fruit - R_fruit);
dR_leaf_dt = ((dC_leaf_dt + MC_leaf_prune)./C_leaf - R_leaf);
dR_stem_dt = (dC_stem_dt./C_stem - R_stem);

% Return all temperature derivatives
dT = cat(1, dT_c_dt, dT_i_dt, dT_v_dt, dT_m_dt, dT_p_dt, dT_f_dt, ...
    dT_s1_dt, dT_s2_dt, dT_s3_dt, dT_s4_dt, ...
    dT_vmean_dt, dT_vsum_dt, ...
    dC_w_dt, dC_c_dt, dC_buf_dt, dC_fruit_dt, dC_leaf_dt, dC_stem_dt, ...
    dR_fruit_dt, dR_leaf_dt, dR_stem_dt);

end