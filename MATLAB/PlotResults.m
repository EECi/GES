%% Plot Results

clear; close all;

parameters;
load out;
load climate;

% Temperature
figure('WindowStyle','docked');
dT=1/(24*60*60/deltaT);
tplot=dT:dT:sim_days;
[~,lenc]=size(tplot);

subplot(3,1,1);
plot(out(:,1)/(24*3600),(out(:,3)-273.15),'b');
hold on;
plot(out(:,1)/(24*3600),(out(:,4)-273.15),'g--');
plot(tplot,climate(1:lenc,1),'k:')
legend on;
legend('Air','Plants','External');
title('Temperatures');


subplot(3,1,2)
plot(out(:,1)/(24*3600),(out(:,2)-273.15),'b');
hold on;
plot(out(:,1)/(24*3600),(out(:,7)-273.15),'c');
plot(out(:,1)/(24*3600),(out(:,8)-273.15),'r');
plot(out(:,1)/(24*3600),(out(:,9)-273.15),'r--');
plot(out(:,1)/(24*3600),(out(:,10)-273.15),'r-.');
plot(out(:,1)/(24*3600),(out(:,11)-273.15),'r:');
plot(tplot,climate(1:lenc,1),'k:')

legend on;
legend('Cover','Floor','Soil (1)','Soil(2)','Soil(3)','Soil(4)','External');
ylabel('Temperature(^oC)');


subplot(3,1,3)
plot(out(:,1)/(24*3600),(out(:,5)-273.15),'m');
hold on;
plot(out(:,1)/(24*3600),(out(:,6)-273.15),'b--');
plot(tplot,climate(1:lenc,1),'k:')
legend on;
legend('Mat','Tray','External');
xlabel('Day');


% Moisture content
figure('WindowStyle','docked');
[len2,~]=size(climate);

for ii=1:len2
C_we(ii)=(climate(ii,4)/100)*sat_conc(climate(ii,1)+273.15);
end

plot(out(:,1)/(24*3600),out(:,14),'c');
hold on;

plot(tplot,C_we(1:lenc),'k:');

legend('Simulation','External');
grid on;
title('Air Moisture Content');
xlabel('Day');
ylabel('Air Moisture Content (kg/m^3)');

% RH
figure('WindowStyle','docked');

for ii=1:len2
RH(ii)=(climate(ii,4)/100);
end

plot(out(:,1)/(24*3600),(out(:,14)./sat_conc(out(:,3))),'c');
hold on;

plot(tplot,RH(1:lenc),'k:');

legend('Simulation','External');
grid on;
title('Relative Humidity');
xlabel('Day');
ylabel('Relative Humidity');


% CO2
figure('WindowStyle','docked');

C_c=out(:,15).*out(:,3)*8.314/(0.044*1.013e5)*1e6;
plot(out(:,1)/(24*3600),C_c,'c');
title('CO_2 concentration');
xlabel('Day');
ylabel('CO_2 concentration (ppm)');

% LAI
figure('WindowStyle','docked');

C_leaf = out(:,18);
C_fruit = out(:,17);
C_stem = out(:,19);

LAI = C_leaf*SLA;

plot(out(:,1)/(24*3600),LAI,'g');
title('Leaf Area Index');
xlabel('Day');
ylabel('Leaf Area Index');