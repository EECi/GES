
% -----------------------------------------------------------------------
%  Greenhouse Energy Simulation model
%  Single zone greenhouse, unheated
% -----------------------------------------------------------------------

clear; close all;

parameters; % Read in all parameters for model from parameters.m

global climate;

% Interpolate weather data

climate=climterp_linear(deltaT); % interpolates climate data

tf = 86400*sim_days;
tspan = [0, tf];

% Initialise data

st = [20 + T_k % Temperature of cover
    12 + T_k % Temperature of air
    12 + T_k % Temperature of vegetation (instantaneous)
    12 + T_k % Temperature of mat
    12 + T_k % Temperature of tray
    12 + T_k % Temperature of floor
    12 + T_k % Temperature of soil Layer 1
    12 + T_k % Temperature of soil Layer 2
    12 + T_k % Temperature of soil Layer 3
    11 + T_k % Temperature of soil Layer 4
    12 + T_k % 24 hour mean temperature of vegetation 
    0 % Temperature sum 
    climate(1,4)/100*sat_conc(climate(1,1) + T_k) % Density of water vapour
    (4e-4)*M_c*atm/(R*(climate(1,1) + T_k)) % Density of CO2
    0.01 % Mass of carbohydrate per unit area of cultivated floor for buffer
    0.001 % Mass of carbohydrate per unit area of cultivated floor for fruit
    0.01 % Mass of carbohydrate per unit area of cultivated floor for leaves
    0.01 % Mass of carbohydrate per unit area of cultivated floor for stem
    0 % Relative growth rate
    0 % Relative growth rate
    0]; % Relative growth rate
    
[t, T] = ode15s(@derivatives, tspan, st);
out = [t, T];
save('out.mat', 'out');

x = 'complete';

