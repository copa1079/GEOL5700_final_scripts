%% 1-D_heat_diffusion_model - MATLAB_R2015b script
%% new code for instanateous step heating of the surface
%% Cole C. Pazar, January 2016, GEOL 5700, Model 2–FTCS

% READ: in order for the code to plot the ATI Alaskan Arctic profiles what
% needs to happen is the data file:   %%%%% (ATI_data_compact.xlsx) %%%%%
% must be imported into matlab by clicking the home tab at the top left and
% importing the data, when importing, use the names Temperature# and Depth#
% respectivley for variable names, to keep the terms consistent with the
% plot for figure(1)

% clear all % un-comment this line if u want the imported data/vars erased

%% constants 
rho = 2000; % kg m^-3 density of the rock
Cp = 2000;
Qm = 0.045; % mantle heat flow
k = 2; % thermal conductivity
kappa = 10^-6; % thermal diffusivity
Ts_i = -20; % initial surface temperature
Ts_f = -7; % final surface temperature

Tstep = Ts_f - Ts_i; % step change in surface temperature

dTdZ = 5/100; % geothermal gradient for the ATI data

%% temporal array

tmin = 0; % initial time in years
tmax = 100; % maximum time in years

dt = 1/365.25; % time step of 1 day

time = tmin:dt:tmax; % time array in years

t = time*3600*24*365.25; % time array in seconds

%% depth array

Zsurface = 0; % surface depth
Zmax = 600; % max depth
dZ = 1; % depth step

Z = Zsurface:dZ:Zmax; % depth array

geotherm = Ts_i + dTdZ.*Z; % linear basic geotherm

T = geotherm+((Ts_f-geotherm).*(1-erf(Z./(2*((kappa*100*3600*24*365.25)^(1/2))))));

% ^ thermal profile for the boundary conditions given in problem

%% reminder: ATI data must be imported now before plotting

% correction variables for the ATI thermal data due to the hot disel in the
% borehole after drilling the hole 

n = 5.093;

c = 6.4-n;
c1 = 6-n;
c2 = 5.6-n;
c3 = 5.5-n;
c4 = 5.4-n;
c5 = 5.21-n;
c6 = 5.096-n;
c7 = 5.093-n;
c8 = 5.093-n;

%% plotting the first figure

figure(1)
    clf

plot(geotherm,Z,'k--','linewidth',2)
    hold on
plot(T,Z,'b','linewidth',2)

% the below set of plots can be removed if you do not have the data file
% containing the ATI profiles, i.e.   ATI_data_compact.xlsx

plot(Temperature-c,Depth,'r--','linewidth',1.5) 
plot(Temperature1-c1,Depth1,'linewidth',1.5)
plot(Temperature2-c2,Depth2,'linewidth',1.5)
plot(Temperature3-c3,Depth3,'linewidth',1.5)
plot(Temperature4-c4,Depth4,'linewidth',1.5)
plot(Temperature5-c5,Depth5,'linewidth',1.5)
plot(Temperature6-c6,Depth6,'linewidth',1.5)
plot(Temperature7-c7,Depth7,'linewidth',1.5)
plot(Temperature8-c8,Depth8,'linewidth',1.5)

set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
axis([-20 10 0 600])
title('Thermal profiles on the north shore of Alaska','fontname','arial','fontsize',18)
ylabel('Depth (m)','fontname','arial','fontsize',18)
xlabel('Temperature (ºC)','fontname','arial','fontsize',18)

%% plotting the rest of the figures

figure(2)
    clf
plot(Temperature,Depth,'r--','linewidth',2)
    hold on
% plot(Temperaturezero,Depthzero,'k--','linewidth',2)
grid on; % sets up  grid

% x0=50; % size
% y0=0; % size
% width=600; % size
% height=800; % size
% set(gcf,'units','points','position',[x0,y0,width,height]) % size

set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
axis([-10.25 6.25 0 500])
title('Arctic subsurface temperature profiles over time: ATI site, AK','fontname', 'arial', 'fontsize', 18)
xlabel('Temperature (ºC)','fontname','arial','fontsize',18)
ylabel('Depth(m)','fontname','arial','fontsize',18)
    pause(3)
plot(Temperature1,Depth1,'linewidth',2)
    pause(1.5)
plot(Temperature2,Depth2,'linewidth',2)
    pause(1.5)
plot(Temperature3,Depth3,'linewidth',2)
    pause(1.5)
plot(Temperature4,Depth4,'linewidth',2)
    pause(1.5)
plot(Temperature5,Depth5,'linewidth',2)
    pause(1.5)
plot(Temperature6,Depth6,'linewidth',2)
    pause(1.5)
plot(Temperature7,Depth7,'linewidth',2)
    pause(1.5)
plot(Temperature8,Depth8,'linewidth',2)
    pause(1.5)
    
    hold off
    
%% corrected figure 

figure(3)
clf

% correction variables

n = 5.093

c = 6.4-n;
c1 = 6-n;
c2 = 5.6-n;
c3 = 5.5-n;
c4 = 5.4-n;
c5 = 5.21-n;
c6 = 5.096-n;
c7 = 5.093-n;
c8 = 5.093-n;

plot(Temperature-c,Depth,'r--','linewidth',2)
hold on
grid on; % sets up  grid

set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
axis([-10.25 6.25 0 500])
title('Arctic subsurface temperature profiles over time: ATI site, AK','fontname', 'arial', 'fontsize', 18)
xlabel('Corrected Temperature (ºC)','fontname','arial','fontsize',18)
ylabel('Depth(m)','fontname','arial','fontsize',18)
    pause(3)
plot(Temperature1-c1,Depth1,'linewidth',2)
    pause(1.5)
plot(Temperature2-c2,Depth2,'linewidth',2)
    pause(1.5)
plot(Temperature3-c3,Depth3,'linewidth',2)
    pause(1.5)
plot(Temperature4-c4,Depth4,'linewidth',2)
    pause(1.5)
plot(Temperature5-c5,Depth5,'linewidth',2)
    pause(1.5)
plot(Temperature6-c6,Depth6,'linewidth',2)
    pause(1.5)
plot(Temperature7-c7,Depth7,'linewidth',2)
    pause(1.5)
plot(Temperature8-c8,Depth8,'linewidth',2)
    pause(1.5)
    
hold off

%% corrected figure close up

figure(4)
clf

plot(Temperature-c,Depth,'r--','linewidth',2)
hold on
grid on; % sets up  grid

set(gca,'fontsize',16,'fontname','arial','YDIR','reverse')
axis([-10.5 -8.5 0 100])
% title('Arctic subsurface temperature profiles over time: ATI site, AK','fontname', 'arial', 'fontsize', 18)
xlabel('Corrected Temperature (ºC)','fontname','arial','fontsize',16)
ylabel('Depth(m)','fontname','arial','fontsize',16)

    pause(3)
plot(Temperature1-c1,Depth1,'linewidth',2)
    pause(1.5)
plot(Temperature2-c2,Depth2,'linewidth',2)
    pause(1.5)
plot(Temperature3-c3,Depth3,'linewidth',2)
    pause(1.5)
plot(Temperature4-c4,Depth4,'linewidth',2)
    pause(1.5)
plot(Temperature5-c5,Depth5,'linewidth',2)
    pause(1.5)
plot(Temperature6-c6,Depth6,'linewidth',2)
    pause(1.5)
plot(Temperature7-c7,Depth7,'linewidth',2)
    pause(1.5)
plot(Temperature8-c8,Depth8,'linewidth',2)
    pause(1.5)
    
    
    hold off
    

%% NOW:    INITIALIZE THE FTCS MODEL

% NEW CONSTANTS:

N = 3; % number of nodes
dz = 10; % 10 meter depth step
dt = 1; % one second time step

% NEW ARRAYS:
T_new = zeros(N,1); 
dTdz = zeros(N,1);
q = zeros(N,1); 

T_new(1) = -10; % surface temperature at node 1
dTdz(N) = 0.025; % heat flux from the mantle [=] ºC m^-1

%% RUN

% calculate T gradient
dTdz(1:N-1) = diff(T_new)/dz;

% calculate heat flux
q = -k*dTdz;

% rate of change of T_new
dqdz = diff(q)/dz;

% update T_new
T_new(2:N) = T_new(2:N) - (1/rho*Cp)*dqdz*dt;

