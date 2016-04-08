%% written by Cole C. Pazar
clear global
clearvars
figure(1)
clf
step = 100;
font = 15;
g = 9.81;

%% parameters

% sand geometry
grain_s = 0.1; % grain size in mm
D = grain_s/1000; % grain diameter in meters

% density
rho_grains = 2650; % density of the grains
rho_air = 1.225; % density of air
porosity = 0.30; % porosity
rho_bulk = rho_grains*(1-porosity)+rho_air*porosity; % bulk density

% particle motion
v_splash = 0.1; % splash velocity of 1/10th of a meter per second
alpha = 10; % angle of impact in degrees
a = alpha*(pi/180); % in radians

% x-array in meters
x_max = 1000;
x_min = 0;
dx = x_max/step; % same size as the grains?
x = (dx/2):dx:x_max-(dx/2); % choose the cell middles

% time step
t_max = 10000; % maximum amount of time
dt = t_max/step; % time step equals the impact number
t = 0:dt:t_max-1;

% plotting controls
imax = length(t);
n_plots = t_max/step;
t_plot = t_max/n_plots;
n_frame = 0;

for i=1:imax;
    
    z = zeros(1,length(x));
    % impact velocity
    z_impact = tand(a)*rand(1,length(x)); % random elevations
    v_impact = 1; % m/s
    % flux law
    dqdx = z_impact*v_impact;
    % topography
    dzdt = -(1/rho_bulk) + dqdx;
    % update
    z = z(i) + dzdt * dt;
        
    if(rem(n_frame,t_plot)==0)
        
    figure(1)
    plot(x,z(i))
    hold on
    plot(x,z)
    axis([0 1000 0 2])
    pause(0.01)
    drawnow 
    
    end
    
    n_frame = n_frame+1;

end

