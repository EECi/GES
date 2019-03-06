function a = sat_conc(T)

T = T - 273.15;

spec_hum = exp(11.56 - 4030./(T + 235)); % Ratio of mass of water vapour to mass of dry air (ref?)
air_dens = -0.0046*T + 1.2978; % Density of air [kg/m^3]
a = spec_hum.*air_dens; % Density of water vapour [kg/m^3]

end