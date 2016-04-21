% go to: https://github.com/copa1079/GEOL5700_final_project   for more outputs etc.
% FINAL PROJECT, model 11
%% UPDATED CODE FOR CORAL GROWTH MODEL         COLE C. PAZAR FOR: GEOL 5700
% eustatic/tectonic coral reef growth
% UPLIFTING TOPOGRAPHY + CHANGING SEA LEVEL
% HUON PINENSULA APPLICATION

% UPDATED: April 20th, 2016.

%% Initialize

    clear global
    clearvars
    figure(1)
    clf
    
    YRs = 3600*24*365; % year in seconds
    step = 150;
    font = 18;
    
    % temporal array
    tmax_yr = 150;             % tmax in ka, 150 ka
    dtyears = 1000;            % time step in years
    dt = YRs*dtyears;          % time step in seconds
    t_max = YRs*tmax_yr*1000;  % tmax in seconds
    t = t_max:-dt:0;           % time array
    
    % spatial array
    x_max = 5000; % max horizontal space in meters
    dx = x_max/step;     % horizontal step in meters
    x = 0:dx:x_max;
    
    % sea level 
    deltaS = 120;    % amplitude of sea level change in meters
    Pyears = 150;     % period in ka
    P = YRs*1000;  % period in seconds
    
    load sealevel_ka.txt
    sea_level = fliplr(transpose(sealevel_ka(:,1)));
    
    % topography
        
    mm = 2.5; % mm for uplift/subsidence rate %%%%%%% can be pos or neg.!!!!
    correction = 2.0;
    D = correction*(mm/1000)*10^-10; % subsidence rate [=] m/s
    m = 0.1;                  % slope of ramp !!!
    topo_max = 180; % maximum height of ramp above mean sea level in meters
    Bmax = -m*x_max+topo_max-(D.*t); % maximum depth of rock
    
    % coral growth
    Gmax = 5; % max growth rate in mm/yr
    Gm = Gmax*3.17*10^(-11); % max growth rate in m/s
    k = 0.1;   % extinction coefficient in m-1
    I_o = 2000; % surface light intensity in microE/m2 s
    I_k = 250;  % saturating light intensity in microE/m2 s
    I = I_o/I_k; % light intensity ratio    
    
    % inital conditions
        N = (x_max/dx)+1; % number of nodes
        
        % sea level
        SL = zeros(N:1); % sea level array
        SL(1:N) = 0;     % initial condition
        
        % topo
        b = -m.*x+topo_max-(D*dt); % inital ramp
        
        % height of coral arrays
        H = zeros(N:1);   % make array for heights of coral relative to SL
        H(1:N) = b;       % initial condition
    
        HD = zeros(N:1);  % make array for height of coral with subsidence
                          % relative to sea level
        HD(1:N)=H-(D*dt); % inital condition
        
        Z = zeros(N:1); % make array for absolute depth of coral below SL
        Z(1:N) = SL-HD; % initial condition
            
        % plotting
        
        imax = length(t);
        
       % n_plots = 150; % how many plots?
       
        n_plots = 3; % how many plots?
        
        t_plot = t_max/n_plots;
     
%% Run

for i = 1:imax

    % sea level with time
    SL(1:N) = sea_level(i); % in meters

    % calculate coral height with subsidence change relative to sea level
    HD = H-(D*t(i)); % in meters
    
    % calculate an absolute depth below current sea level
    Z = SL-HD; % in meters
    
    % calculate the growth of the coral for this time step
    dh = (Gm*tanh(I*exp(-k*Z)))*dt; % in meters
    
    % correct dh for negative Z values so there is no growth
    neg = find(Z<0); % find spots where coral is above sea level
    dh(neg)=0; % corrects dh array with 0 growth for coral above sea level
         
    % calculate the new height of the coral
    H(1:N) = HD + dh; % in meters
    
    % bedrock topography
    uplift = (D*t(i)); % try to add it in
    B = -m.*x+topo_max-uplift; % in meters
     
    % plotting
    if  (rem(t(i),t_plot)==0)

    figure(1)
    clf
    hold all    
    fill(x/1000,smooth(smooth((H))),'g','linewidth',1.5) % corals !!!
    plot(x/1000,SL,'b','linewidth',2)
    plot(x/1000,zeros(length(x))+9,'r','linewidth',1.5) % max sea level
    plot(x/1000,zeros(length(x))-122,'r','linewidth',1.5) % min sea level
    xlabel('distance along profile (km)', 'fontname', 'arial',...
        'fontsize',font)
    ylabel('elevation relative to zero in sea level (m)',...
        'fontname', 'arial', 'fontsize',font)
    title('Coral reef growth model, with actual scaled global mean sea level',...
        'fontname', 'arial', 'fontsize',font)
    set(gca,'fontsize',font,'fontname','arial')
    ht=text(4,-200,['time: ',num2str((t(i)/YRs/1000)), ' ka'],...
        'fontsize',font); % print time in animation
    ht2=text(4,-220,['sea level: ',num2str(round(SL(i))), ' m'],'fontsize',font);
    ht3=text(4,-240,[num2str(mm),' mm/yr subsidence rate'],'fontsize',font); 
    ht4=text(4,-260,[num2str(Gmax),' mm/yr max growth rate'],'fontsize',font); 
    ht5=text(4,-280,[num2str(m),' m/km slope'],'fontsize',font); 
    legend('corals','sea level','sea level range')
    axis([0 (x_max/1000) -400 100])
    hold all
    
    pause(0.01)
    
    end
end

    % spatial array
    x_maxx = 150; % max horizontal space in meters
    dxx = x_maxx/step;     % horizontal step in meters
    xx = x_maxx:-dxx:0;

    figure(2)
    clf
    plot(xx,sea_level,'ko','linewidth',4)
    hold on
    plot(xx,sea_level,'b','linewidth',1.5)
    plot(xx,smooth(sea_level),'c','linewidth',2)
    plot(xx,ones(length(xx)),'k','linewidth',1)
    xlabel('time (ka)', 'fontname', 'arial','fontsize',font)
    ylabel('mean sea level realtive to today (m)','fontname','arial',...
        'fontsize',font)
    title('Scaled mean global sea level curve, from deep sea \delta18O'...
    ,'fontname', 'arial', 'fontsize',font)
    set(gca,'fontsize',font,'fontname','arial')
    legend('\delta18O points','sea level','smoothed curve','modern sea level')
    axis([0 150 -140 60])

