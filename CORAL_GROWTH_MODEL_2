%% UPDATED CODE FOR CORAL GROWTH MODEL
%% UPLIFTING TOPOGRAPHY + CHANGING SEA LEVEL
%% COLE C. PAZAR GEOL 5700
% load Pleist_oxy.txt
% load pleist_del18O.txt
% load ox_200.txt
% load yrs_200.txt
%% Initialize
clear all
YRs = 3600*24*365; % year in seconds

    % temporal array
    dtyears = 200;          % time step in years
    dt = YRs*dtyears;       % time step in seconds
    tmax_yr = 250000;       % tmax in years
    tmax = YRs*tmax_yr;     % tmax in seconds
    t = 0:dt:tmax;          % time array
    
    % spatial array
    dx = 10;     % horizontal step in meters
    xmax = 3000; % max horizontal space in meters
    x = 0:dx:xmax;
    
    % sea level 
    deltaS = 60;     % half amplitude of sea level change in meters
    Pyears = 120000; % 120 ka period in years
    P = YRs*Pyears;  % period in seconds
    
    % Topography
   
    mm = -2; % mm for uplift/subsidence rate %%%%%%% can be pos or neg.!!!!
    
    D = (mm/1000)/(P/Pyears); % subsidence rate [=] m/s
    m = 0.5;                  % slope of ramp !!!
    topo_max = 200; % maximum height of ramp above mean sea level in meters
    Bmax = -m*xmax+topo_max-(D*tmax); % maximum depth of rock
    
    % Coral growth
    
    Gmax = 10; % max growth rate in mm/yr
    Gm = Gmax*3.17*10^(-11); % max growth rate in m/s
    k = 0.1;   % extinction coefficient in m-1
    Io = 2000; % surface light intensity in microE/m2 s
    Ik = 250;  % saturating light intensity in microE/m2 s
    I = Io/Ik; % light intensity ratio
    
    % Inital Conditions
        N = (xmax/dx)+1; % number of nodes
        
        % sea level
        
        SL = zeros(N:1); % sea level array
        SL(1:N) = 0;     % initial condition
        
        % topo
        b = -m.*x+topo_max; % inital ramp
        
        % Height of coral arrays
        
        H = zeros(N:1); % make array for heights of coral relative to SL
        H(1:N)=b;       % initial condition
    
        HD = zeros(N:1);  % make array for height of coral with subsidence
                          % relative to sea level
        HD(1:N)=H-(D*dt); % inital condition
        
        Z = zeros(N:1); % make array for absolute depth of coral below SL
        Z(1:N) = SL-HD; % initial condition
 
    imax = length(t/10);
    
%% Run

for i=1:1:imax

    % sea level with time
    SL(1:N) = deltaS*sin((2*pi*t(i))/P); % in meters
    
    % Calculate coral height with subsidence change relative to sea level
    HD = H-(D*dt); % in meters
    
    % Calculate an absolute depth below current sea level
    Z = SL-HD; % in meters
    
    % Calculate the growth of the coral for this time step
    dh = (Gm*tanh(I*exp(-k*Z)))*dt; % in meters
    
    % correct dh for negative Z values so there is no growth
    neg = find(Z<0); % find spots where coral is above sea level
    dh(neg)=0; % corrects dh array with 0 growth for coral above sea level
         
    % Calculate the new height of the coral
    H(1:N) = HD + dh;
    
    % bedrock topography (linear for now)
    B = -m.*x+topo_max-D*t(i); % in meters
     
    % plotting
    figure(1)
    clf
    hold all
    plot(x,B,'k','linewidth',2)  
    plot(x,SL,'b','linewidth',1.5)
    plot(x,zeros(length(x))+60,'r--','linewidth',1.5)
    plot(x,zeros(length(x))-60,'c--','linewidth',1.5)
    fill(x,H,'g')
    grid on
    xlabel('distance (m)', 'fontname', 'arial', 'fontsize',18)
    ylabel('depth (m)', 'fontname', 'arial', 'fontsize',18)
    title('Coral reef growth model','fontname', 'arial', 'fontsize',18)
    set(gca,'fontsize',16,'fontname','arial')
    ht=text(1050,350,['time: ',num2str(round(t(i)/YRs/1000)), ' ka  '],...
        'fontsize',16); % print time in animation
    ht2=text(1050,150,['sea level change amplitude = 120 m'],'fontsize',16); 
    ht3=text(1050,250,['2 mm/yr uplift rate'],'fontsize',16); 
    legend('topography','sea level')
    axis([0 1500 -400 400])
    hold all
    pause(0.1)
    
end

