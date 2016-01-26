%% THERMAL PROFILES WRITTEN JAN 13TH, 2016 IN GEOL5700 - COLE C. PAZAR
%  ASSUMING NO RADIOGENIC HEAT COMPONENT
clearvars
clear all
figure(1)
clf

%% CONSTANTS

Ts_bar = -10; % mean annual surface temperature 
k = 2.5;
Qm = 0.045;
kappa = 10^-6;
dT = 15; % amount of surface oscillation [=] ºC

P_day = 24*60*60; % daily period of oscillation [=] seconds
P_annual = P_day*365; % annual period of oscillation [=] seconds


%% ARRAYS

dz = 1;
zmax = 15;
z = 0:dz:zmax; % depth array [=] meters
zmid = z(1:end-1)+(dz/2);

dt = 5*60*60*24;% 5 days [=] seconds
tmax = 1*P_annual;
time = 0:dt:tmax;
imax = length(time);

zstar_annual = sqrt(kappa*P_annual/pi); % approximate annual e-folding depth [=] meters

Tzero = zeros(size(z)); 

z1 = 0.05;
z2 = 0.11;
z3 = 0.31;

z_1 = find(zmid==z1);
z_2 = find(zmid==z2);
z_3 = find(zmid==z3);

%%  PART 2) ANNUAL OSCILLATIONS

% RUN

for i = 1:imax
    T = Ts_bar + dT*exp(-z./zstar_annual).*sin((2*pi*time(i)/P_annual)-(z./zstar_annual));
    
    Tmax = Ts_bar + dT*exp(-z./zstar_annual);

    Tmin = Ts_bar - dT*exp(-z./zstar_annual);
    
figure(1)
plot(Tzero,z,'c--','linewidth',2)
title('Annual-scale oscillation of subsurface temperatures', 'fontname', 'arial', 'fontsize', 18)
xlabel('Temperature (ºC)','fontname','arial','fontsize',18)
    ylabel('Depth below the surface (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
    axis([Ts_bar-dT-1 Ts_bar+dT+1 0 zmax])
hold on
plot(Tmax,z,'r','linewidth',2)
plot(Tmin,z,'r','linewidth',2)
pause(0.1)
plot(T,z,'k')


end
