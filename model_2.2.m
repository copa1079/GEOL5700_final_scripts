%% Borehole Thermal Profiles in Cape Thompson, Alaska
% Incorporates the data from alaska, best in matlab versions from 2015+
% Finite Difference Model FTCS Diffusion
% Cole C. Pazar, GEOL5700, January 2016

% parameter optimization done with the chi-squared value

%% modeling the rest of the data
clear all

figure (1) % thermal model with data
clf
figure(2) % chi-squared
clf

% requires 8 data files:
% download the below data in a zip file, must be in the same folder as this
                                                       % matlab code script
load AWU_81AUG22.txt
load AWU_82SEP06.txt
load AWU_83AUG30.txt
load AWU_84JUL29.txt
load AWU_89JUL31.txt
load AWU_02AUG14b.txt
load AWU_08JUL24.txt
load AWU_12AUG07.txt

z81 = AWU_81AUG22(:,1); % grabs the data from column 1
T81 = AWU_81AUG22(:,2); % grabs the data from column 2

z82 = AWU_82SEP06(:,1); 
T82 = AWU_82SEP06(:,2); 

z83 = AWU_83AUG30(:,1); 
T83 = AWU_83AUG30(:,2);

z84 = AWU_84JUL29(:,1); 
T84 = AWU_84JUL29(:,2); 

z89 = AWU_89JUL31(:,1); 
T89 = AWU_89JUL31(:,2); 

z02 = AWU_02AUG14b(:,1);
T02 = AWU_02AUG14b(:,2);

z08 = AWU_08JUL24(:,1);
T08 = AWU_08JUL24(:,2);

z12 = AWU_12AUG07(:,1); 
T12 = AWU_12AUG07(:,2); 

numdata12 = length(T12); % number of data points in column 2 for 2012

% dTdz_base = mean(diff(Tdata(end-10:end))./diff(zdata(end-10:end))); 
% ^ calculates the geothermal gradient at the base of the data

% material properties

Qm = 0.05;
k = 1.6; % W/(m*K) thermal conductivity
rho = 2700; % density of the rock
c = 2184; % heat capacity
kappa = k/(rho*c);

dTdz_base = Qm/k; % ºC/km

% depth array
dz = 1;
zmax = 400;
z = 0:dz:zmax;
N = length(z);

% time arrays and controls
ndays = 365*150;
dt = 24*3600; % in seconds; again take charge of the important time step
tmax = ndays*3600*24; % max time in Seconds
t=0:dt:tmax; % time array: note this has nothing to do with N which is a 
             % number of spatial array elements
imax = length(t); % note imax and N are completely independent

% plotting controls

nplots = 100;
tplot = tmax/nplots;
 
dTdz0 = dTdz_base; % sets bottom bc gradient

Ts_old = -9; % initial top boundary condition
Ts_new = -1; % final boundary top condition     % ?T = 8º C or 8K
T = Ts_old + (dTdz0*z); % this is the whole profile at t=0
T0=T; % saves initial profile for plotting
q = zeros(1,N); % fills an array the same size as T wiht zeros.
%in reality this is shifted downward by 1/2 node:fluxes in bottoms of boxes
j=0; % size of chi-squared values ?

%% run loop

for i=1:imax
    
    T(1) = Ts_new; % this way you can take charge of this wiht any
    % prescribed history of Ts over time... like a sinsoid, etc...

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

% correction variables for the AWU data set:

Ti = 3.33; % temperature at 400 meters for all profiles, set to 2002 data

c81 = 6.698-Ti;
c82 = 5.554-Ti;
c83 = 4.090-Ti;
c84 = 3.902-Ti;
c89 = 3.477-Ti;
c02 = 3.330-Ti;
c08 = 3.330-Ti;
c12 = 3.330-Ti;

% plotting

if (rem(t(i),tplot)==0)
    j = j+1;
    timeplot(j) = t(i);
    
    figure(1)
 plot(T,z,'k','linewidth',2)
    hold on
 plot(T0,z,'k-','linewidth',1)
 
 plot(T81-c81,z81,'linewidth',1)
 plot(T82-c82,z82,'r','linewidth',1)
 plot(T83-c83,z83,'y','linewidth',1)
 plot(T84-c84,z84,'m','linewidth',1)
 plot(T89-c89,z89,'c','linewidth',1)
 plot(T02-c02,z02,'y.','linewidth',1)
 plot(T08-c08,z08,'m.','linewidth',1)
 plot(T12-c12,z12,'c.','linewidth',1)

 legend('thermal model','linear geotherm','initial data, 1981'...
     ,'data 1982','data 1983','data 1984','data 1989','data 2002'...
     ,'data 2008','data 2012','Location','northeast') % creates a legend
        grid on
        set(gca,'YDIR','reverse')
 title('Northern AK permafrost temperatures: AWU site','fontname',...
     'arial', 'fontsize', 16)
    axis([-10 5 0 300])
        xlabel('Temperature [corrected] (°C)','fontname','arial','fontsize',16)
        ylabel('Depth below the surface (m)','fontname','arial','fontsize',16)
        set(gca,'fontsize',16,'fontname','arial')
    drawnow
 hold off
  
 Tmodel = interp1(z,T,z12);  % picks the points in the model at a specific 
                             % depth, allowing for same size arrays
                              
 chisq(j) = (1/numdata12)*sum((T12-Tmodel).^2);  % statistical check
                                                 % missing the 1/sigma **
  
 % minimum value near ~98 years since an instantaneous thermal
 % step change heating problem

end

end

%% finalize plots

figure(2)
plot(timeplot/(3600*24*365),chisq,'linewidth',2)
title('Statistical fit: transient step-heating to -1 ºC','fontname',...
    'arial', 'fontsize', 16)
    axis([60 150 0 0.1])
    grid on
        xlabel('Time (years)','fontname','arial','fontsize',16)
        ylabel('X^2 (chi-squared value)','fontname','arial','fontsize',16)
        set(gca,'fontsize',16,'fontname','arial')

%% end of code
