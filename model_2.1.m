%% Borehole Thermal Profiles in Cape Thompson, Alaska
% This code takes real data and fits the thermal perturbation to the 
% borehole data (figure 1), doing a chi-squared regression (figure 2)
% Finite Difference Model FTCS Diffusion
% Cole C. Pazar, GEOL5700, January 2016
% with help from Robert S. Anderson

%% Open Figures
clear all
figure (1) % thermal model with data
clf
figure(2) % chi-squared
clf

load cape_thompson.dat % loads the arctic data, must be in the same file

zdata = cape_thompson(:,1); % grabs the data from column 1
Tdata = cape_thompson(:,2); % grabs the data from column 2
numdata = length(Tdata); % number of data points in column 2

% dTdz_base = mean(diff(Tdata(end-10:end))./diff(zdata(end-10:end))); 
% ^ calculates the geothermal gradient at the base of the data

dTdz_base = 0.02; % ºC/km

% depth array
dz = 1;
zmax = 400;
z = 0:dz:zmax;
N = length(z);

% material properties
k = 2; % W/(m*K) thermal conductivity
rho = 2000; % density of the permafrost
c = 2000; % heat capacity
kappa = k/(rho*c);

% time arrays and controls
ndays = 365*150;
dt = 24*3600; % in seconds; again take charge of the important time step
tmax = ndays*3600*24; % max time in Seconds
% dt = tmax/N; % time step ?
t=0:dt:tmax; % time array: note this has nothing to do with N which is a 
             % number of spatial array elements
imax = length(t); % note imax and N are completely independent

% Plotting controls

nplots = 100;
tplot = tmax/nplots;

% dTdz0 = 0.025; % sets bottom bc gradient
dTdz0 = dTdz_base;

Ts_old = -7.1; % initial top boundary condition
Ts_new = -5.2; % final boundary top condition
T = Ts_old + (dTdz0*z); % this is the whole profile at t=0
T0=T; % saves initial profile for plotting
q = zeros(1,N); % fills an array the same size as T wiht zeros.
%in reality this is shifted downward by 1/2 node:fluxes in bottoms of boxes
j=0; % size of chi-squared values ?

%% run loop

for i=1:imax
    
    T(1) = Ts_new; % this way you can take charge of this wiht any
    % prescribed history of Ts over time... like a sinsoid, or a
    % ramp...or...

% calculate T gradient
% dTdz(1:N-1) = diff(T)/dz;
dTdz = diff(T)/dz;

% calculate Heat Flux
q(1:end-1) = -k*dTdz;
q(end) = -k*dTdz0; % bottom boundary condition

% rate of change of T
dqdz = diff(q)/dz;

% update T
T(2:N) = T(2:N) - (1/(rho*c))*dqdz*dt;

% plotting

if (rem(t(i),tplot)==0)
    j = j+1;
    timeplot(j) = t(i);
    
    figure(1)
 plot(T,z,'r','linewidth',2)
 hold on
 plot(T0,z,'k--','linewidth',2)
 plot(Tdata,zdata,'ko')
 set(gca,'YDIR','reverse')
 title('Cape Thompson permafrost thermal evolution','fontname','arial',...
     'fontsize',16)
        xlabel('T (°C)','fontsize',16,'fontname','arial')
        ylabel('Depth (m)','fontname','arial','fontsize',16)
        set(gca,'fontsize',16,'fontname','arial')
 drawnow
 hold off
  
 Tmodel = interp1(z,T,zdata); % picks the points in the model at a specific 
                              % depth, allowing for same size arrays
                              
 chisq(j) = (1/numdata)*sum((Tdata-Tmodel).^2);  % statistical check
 
 % minimum value near ~111 years since an instantaneous thermal
 % step change heating problem

end

end

%% finalize plots

figure(2)
plot(timeplot/(3600*24*365),chisq)
title('Statistical fit: min at 111 years','fontname','arial',...
    'fontsize',16)
        xlabel('Time (years)','fontname','arial','fontsize',16)
        ylabel('chi-square value','fontname','arial','fontsize',16)
        set(gca,'fontsize',16,'fontname','arial')
