%% Hydrologic overland flow model, re-done over spring break.
%  Finalized on March 26th, 2016 while flying over Greenland.
%  Cole Pazar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear global
clearvars
figure(1)
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code organization constants
font     = 18;   % choose the font size for plots
line_wth = 3;    % choose the line width
border   = 2;    % above and below the max and min for topography
border_2 = 4;    % choose the border for the rainfall plot
factor   = 4;    % factor for the time step size
step     = 200;  % holding the array lengths to be 200 elements long
day      = 24;   % number of hours in a day
sec      = 60;   % number of seconds in a min
min      = 60;   % number of min in an hr   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize
% physical constants
g = 9.81;  % gravitation
k = 0.4;   % von karman's constant from Anderson and Anderson (2010)
b = 10^-4; % roughness of the bed, 'z_o' in Anderson and Anderson (2010)
           % also known as Nikuradse's roughness. also:  b ~ D84/10
           % this value affects the numerical stability of the code
% set up the x-array
x_min = 0;
x_max = 1000;
dx    = x_max/step;
x     = x_min:dx:x_max; % distance array
% set up the time array
t_min        = 0;
hrs          = day/2;                        % half of a day
t_max        = hrs*min*sec;               % maximum alloted time
time_stepper = factor*1000;               % for the time step
dt           = t_max/(step*time_stepper); % time step
           t = t_min:dt:t_max;            % time array
% randomization
           xx = 10^-3;
random_number = xx*randn(size(t));

% set up a linear topography
S     = -0.01;       % constnat slope of the channel bed
z_min = 0;           % minimum elevation
z_max = 10;          % maximum elevation
z = (S * x) + z_max; % elevation array
% set up 2 rainstorms
initial      = 0.01;                  % initial amount of rain
oscillations = 1.8;                   % choose the number of storms
% minimum at 6 and at the edges
period       = t_max/oscillations;    % choose the period of the storms
amp          = 1;                     % choose the size of the storms
amplitude    = initial*amp;           % actual amplitude
% set up precipitation, P, while updating its values
shift            = -1;
P                = initial+amplitude*sin(2*pi*(t+(shift*(sec*min)))/period)+random_number;
    P            = P/(sec*min); % converting
         P       = max(P,0);
negative         = find(P<0);
P(negative) = 0; % can't have negative amounts of precipitation
% set up infiltration rate
infil = 0.01; % infiltration [=] m/hr
infil = infil/(sec*min); %   [=] m/s
% set up water height arrays
H = zeros(size(x)); % fixes the height variables throughout my code
Z = z+H; % total topography including water thickness
% plotting controls
nplots = day*4;           % number of plots
tplot  = t_max/nplots; % plotting interval
imax   = length(t);    % should be step length
tic % time counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run

for i = 1:imax
 
    % thickness variables to live at the edges, known as the variable 'z' 
    % in Anderson and Anderson (2010)
    H_edge = (H(1:end-1))+(H(2:end));
    
    % set up the infiltration over the x-array
    I = zeros(size(x));
    
    % this finds the wet part of the slope
    water = find(H>0); 
    
    % infiltration rate where there is water
    I(water) = infil;
    
    % uniform rain over the x-array
    R = P(i)*ones(size(x));
        
    % calulates the actual topogaphic slope
    dzdx = diff(z)/dx; 
    
    % LAW OF THE WALL FOR FLUID FLOW
    % set up u_star for cell edges
    U_star = sqrt(g*H_edge.*abs(dzdx)); % square root of shear stress over the fluid density
    % for stability, w finds thicknesses greater than the bed roughness:
    w = find(H_edge>b);
    % set up the average velocity
    U_bar = zeros(size(H_edge));
    U_bar(w) = (U_star(w)./k).*(log(H_edge(w)/b)-1);
    
    % set up the discharge
       Q = U_bar.*H_edge;
       Q = [0 Q];            % fixes the upper B.C.
    dQdx = diff(Q)/dx;       
    dQdx = [dQdx dQdx(end)]; % fixes additional B.C.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conservation statement:
    %
    % time rate of change of water thickness
    % equals the fluxes downslope plus the
    % precipitation rate of rain, minus the
    % infiltraion rate into the ground.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dHdt = -dQdx + ( R - I );% conservation statement
       H = H + (dHdt*dt);    % updates the thickness as the rate varies
       H = max(0,H);         % negative thickness of water is impossible
       Z = z + H;            % updates the total thickness with topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up plots
    if(rem(t(i),tplot)==0) 
    
    figure (1) % animation  
    % thickness plot
    subplot(2,2,1) 
    plot(x,H*100,'b','linewidth',line_wth)
    axis([0 x_max -0.1 1])
    ht=text(800,0.9,['  ',num2str(t(i)/(min)),' min'],'fontsize',font);
    grid on
    set(gca,'fontsize',font,'fontname','arial')  
    title('water depth vs. distance downslope','fontsize',font,'fontname','arial')
    ylabel('flow thickness [cm]','fontsize',font,'fontname','arial')
    xlabel('distance [m]','fontsize',font,'fontname','arial')
    % hillslope plot
    subplot(2,2,2)
    plot(x,z,'k','linewidth',line_wth) 
    axis([0 x_max z_min-border z_max+border])
    grid on
    set(gca,'fontsize',font,'fontname','arial')  
    title('linear hillslope over 1 kilometer','fontsize',font,'fontname','arial')
    ylabel('relative elevation [m]','fontsize',font,'fontname','arial')
    xlabel('distance [m]','fontsize',font,'fontname','arial')
    % flux plot
    subplot(2,2,3) 
    plot(x,Q*1000,'r','linewidth',line_wth)
    axis([0 x_max -0.5 4])
    ht=text(800,3.5,['  ',num2str(t(i)/(min)),' min'],'fontsize',font);
    grid on
    set(gca,'fontsize',font,'fontname','arial')  
    title('flux of water vs. distance downslope','fontsize',font,'fontname','arial')
    ylabel('discharge [mm^2/s]','fontsize',font,'fontname','arial')
    xlabel('distance [m]','fontsize',font,'fontname','arial')
    % storm plot
    subplot(2,2,4) 
    plot(t/(sec*min),P*(sec*min)*1000,'b','linewidth',0.5)
    axis([0 t_max/(sec*min) -border_2 (max(P)*sec*min*1000)+border_2])
    grid on
    set(gca,'fontsize',font,'fontname','arial')  
    title('prescribed storm over 720 mins','fontsize',font,'fontname','arial')
    xlabel('time [hours]','fontsize',font,'fontname','arial')
    ylabel('rain [mm/hr]','fontsize',font,'fontname','arial') 
    drawnow
    pause (0.01)
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
