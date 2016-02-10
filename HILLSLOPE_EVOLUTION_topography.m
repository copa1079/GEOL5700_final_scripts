%% CODE FOR HILLSLOPE EVOLUTION WITH A REAL TOPOGRAPHIC PROFILE, evolving over 1 million years.
%% WRITTEN BY COLE C. PAZAR, ON FEBURARY 8TH, 2016.
%  1D diffusion staggered variables code style inspired by RSA

% now with a real topographic profile
% requires the file named: profile_for_erosive_diffusion.txt to be in the same folder!
% the topographic data is included in the git hub repo.

%  clear all

clearvars
figure(1)
clf

%% initialize

k = 0.1;     %  relatively high hillslope diffusivity [=] m^2/year

% distance array [=] m
step = 1000;
xmin = 0;
xmax = 27844;
dx = xmax/step; % was 1
x = xmin:dx:xmax;

% initial topography and soil thickness [=] m
disp = 0; % amount of topographic displacement
slopedegrees = -25;      % slope 
m = tand(slopedegrees);  % slope re-calculated
zmid = 4000;              % elevation at middle of an already eroded hillslope
% zb = (m.*x) + zmid;      % topography
% zb = -((0.1/500).*(x+disp).^2) + zmid;  % topography

load profile_for_erosive_diffusion.txt
topo = transpose(profile_for_erosive_diffusion(:,1));
zb = transpose(profile_for_erosive_diffusion(:,1));

H = zeros(1,length(x));  % soil thickness
z = zb+H;                % net

% weathering and transport parameters

wnot = 10^-3; % [=] m/yr
wstar = 1;  % scale for weathering rate, 1 meter
A = 0.01;
hstar = 0.2;
wdot = zeros(1,length(x)); % sets up array for weathering rate

% channel incision

edot = 10^-2; % [=] m/yr

% temporal array

tmax = 1000000; % 1 million years
dt = 10;
t = 0:dt:tmax;
nplots = 100;
tplot = tmax/nplots;

% fault w/random variables

dipdegrees = -70;
dip = tand(dipdegrees);
fault = (dip.*(x+300))+zmid-2000;
horzslip = 200;              % on average 20 meters of slip every 1000 years
vertslip = horzslip*(tand(dipdegrees)-tand(slopedegrees)); % accounts for geometry of the fault
slipperiod = 20000;          % fault period: every 10 ka


%% run

for i=1:length(t)

    wdot=wnot.*exp(-H./wstar); % exponential weathering rate
    
       % or :
       % humped weathering rate
       % wdot=wnot.*exp(H./wstar)+((A*H./wstar).*exp(-H./wstar));

    slope = diff(z)/dx;
    
    % find peaks and troughs
    
    xedge = x(1:end-1)+(dx/2);
    lt = find(slope<0); % slopes with left peaks
    rt = find(slope>0); % slopes with right peaks

    hedge = zeros(size(slope)); % sets up matrix sizes for edge variables
    hedge(1) = H(2);
    hedge(end)= H(end-1);
    hedge(rt) = H(rt+1);
    hedge(lt) = H(lt);
    
    q = -k.*slope.*(1-(exp(-hedge./hstar))); % discharge 
    
    q = [q(1) q q(end)]; % boundary conditions
    
    % develop topography and conserve mass
    dh_dt = wdot-diff(q)/dx;
    H = H + (dh_dt*dt);
    H = max(0,H); % for another boundary setup...
    zb = zb - (wdot*dt);
    z = H+zb;
   
    if(rem(t(i),tplot)==0)
        
    figure(1)
    %subplot(2,1,1)
    subplot('position',[0.08 0.4 0.85 0.55])
    plot(x/1000,z,'g','linewidth',3)
    hold on
    plot(x/1000,zb,'k','linewidth',2)
    hold on
    plot(x/1000,topo,'r','linewidth',1.25)
    ht=text(12,4000,['  ',num2str(t(i)/1000), ' ka '],'fontsize',18); % this makes the years print on screen
    axis([0 28 2900 4450]) 
    title('Landscape evolution of a topographic transect, Sawatch Range, Colorado, kappa = 0.1 m^2/year, over 1 Ma','fontsize',20)
    ylabel('relative elevation [m]','fontsize',18)
    set(gca,'fontsize',18)
    
    legend('Regolith','Bedrock','Initial profile','Location','southeast')
    
    hold off 
    %subplot('position',[0.08 0.08 0.85 0.25])
    %RI = imref2d(size(I));
    %RI.XWorldLimits = [0 3];
    %RI.YWorldLimits = [2 5];
    %subplot('position',[0.08 0.4 0.4 0.4])
    
    %subplot(2,1,2)
    subplot('position',[0.08 0.08 0.85 0.25])
    plot(x/1000,H,'g','linewidth',3)
    axis([0 28 0 150])
    xlabel('north-to-south distance along profile[km]','fontsize',18)
    ylabel('regolith thickness [m]','fontsize',18)
    set(gca,'fontsize',18)
    grid on
        
    drawnow
    
    end
   
end

%   this part of the code is in case you'd like to read in
%   the image of the transect across a DEM map of the surface. 

%    figure(2)
%    clf
%    I = imread('topographic_transect_sawatch.jpg');
%    imshow(I);


%% end of code
