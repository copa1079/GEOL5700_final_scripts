%% THERMAL PROFILES WRITTEN JAN 13TH, 2016 IN GEOL5700 - COLE C. PAZAR
%  ASSUMING NO RADIOGENIC HEAT COMPONENT
% clearvars
% clear all
figure(1)
clf

%% CONSTANTS

Ts_bar = -10; % mean annual surface temperature 
k = 2.5; % % W m^-1 ºC^-1
Qm = 0.045; % W m^-2
kappa = 10^-6; % mm^2 s^-1
dT = 15; % amount of annual surface oscillations [=] ºC
dT_daily = 10; % amount of daily surface oscillations [=] ºC

P_day = 24*60*60; % daily period of oscillation [=] seconds
P_annual = P_day*365; % annual period of oscillation [=] seconds


%% ARRAYS

dz = 1;
zmax = 15;
z = 0:dz:zmax; % depth array [=] meters

dz2 = 1/15;
zmax2 = 1;
z2 = 0:dz2:zmax2; % depth array [=] meters


dt = 5*60*60*24;% 5 days [=] seconds
tmax = 1*P_annual;
time = 0:dt:tmax;
imax = length(time);

zstar_annual = sqrt(kappa*P_annual/pi); % approximate annual e-depth [=] meters
zstar_daily = sqrt(kappa*P_day/pi); % approximate annual e-depth

Tzero = zeros(size(z)); 

%% PART 3) NOW FOR OSCILLATORY DAILY-ANNUAL VIDEO

% RUN

for i = 1:imax
    
    T = Ts_bar + dT*exp(-z./zstar_annual).*sin((2*pi*time(i)/P_annual)-(z./zstar_annual)) 
    
    T_full = T + dT_daily*exp(-z2./zstar_daily).*sin((2*pi*time(i)/P_day)-(z2./zstar_daily));

figure(1)

plot(T,z,'k','linewidth',2)
title('subsurface annual thermal range to 15 meters', 'fontname', 'arial', 'fontsize', 18)
xlabel('temperature (ºC)','fontname','arial','fontsize',18)
    ylabel('depth below the surface (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
    axis([Ts_bar-dT-1 Ts_bar+dT+1 0 15])
hold on
plot(T,z,'k')
hold on
plot(Tzero,z,'k--')
hold off

    pause(0.05) % pauses the interval between frame pastes
    
end

