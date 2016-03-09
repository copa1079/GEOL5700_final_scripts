%% OVERLAND FLOW MODEl
%  written by Cole C. Pazar in SI units
clear global
clearvars
figure(1)
clf
% conservation statement:
% change in water height is a function of fluxes, rainfall, infiltration
 
% can update topography as QGIS data, make sure step is the same size!!
 
%% Initialize
 
    step = 1000;
 
% topographic constants
 
    Zmax = 500;   % max elevation of valley in meters
    m    = 0.1;   % slope of valley (on linear topography)
    
% overland flow constants
 
    n = 0.1;
    
% rainfall constants ( for changing the storm as well )
 
    minutes        = 15;                 % duration of the storm in minutes
    stormlength    = 60*minutes;         % length of storm in seconds
    stormvelocity  = 0.5;                % velocity storm moves down slope
    stormrainheight = 0.02/stormlength;
 
% distance array
 
    xmin = 0;         % min horizontal space in meters
    xmax = 2000;      % max horizontal space in meters
    dx = xmax/step;   % horizontal step in meters
    x = xmin:dx:xmax;
    xedge = xmin+(dx/2):dx:xmax-(dx/2);
    
% time array
 
    maxminutes = 200;      % number of minutes to run the simulation
    tmax = 60*maxminutes;  % tmax in seconds
    dt = tmax/step;        % time step in seconds
    t = 0:dt:tmax;         % time array
    
% plotting controls
 
    N = length(x);         % length of all nodes
    
    imax = length(t);      % loop length
    nplot = 2*maxminutes;  % number of plots
    tplot = tmax/nplot;    % times at which plots will be made
    xx = [x, fliplr(x)];
    bottomline=-200000000.*ones(1,N);
    
% ground arrays
 
    earthelevation = zeros(1,N);       % for elevation profile
    earthelevation(1:N) = Zmax - m.*x; % linear valley profile
    slope = zeros(1,N);                % preallocate slope array
    
% water arrays
 
    Hinitial = 0;            % initial water height
    H = Hinitial.*ones(1,N); % create array of water height
    Hedge = zeros(1,N-1);    % create array for height of water at edges
    jmax = length(Hedge);    % loop for filling Hedge array
    for j=1:jmax
        Hedge(j) = (H(j)+H(j+1))/2; % calculate thickness of water at edges
    end 
    waterelevation = earthelevation(1:N) + H(1:N);
    
% flux and strom arrays
    
    Q    = zeros(1,N);           % preallocate Q array
    dQdx = zeros(1,N-2);         % preallocate dQdx array
    
% hydrograph tracking
 
    Hydro_1 = zeros(1,imax);
    Hydro_2 = zeros(1,imax);
    Hydro_3 = zeros(1,imax);
    
% storm arrays
 
    stormfront = zeros(1,imax);
    stormtail  = zeros(1,imax);
   
    
%% Run
 
for i = 2:imax
    
    % slope determination
    
    slope = abs(diff(waterelevation./dx));
    
    % flux from flow
    
    Q(2:N) = 1/n .* Hedge.^(5/3).*slope.^(1/2); % manning equation, average
    Q(1) = 0;                                   % fixes B.C.
    
    Hydro_1(i) = Q(N);
    Hydro_2(i) = Q(ceil(N/2));
    Hydro_3(i) = Q(ceil(N/10));
    dQdx(1:N-1) = diff(Q)./dx; % diff Q to get dQdx
    dQdx(N) = dQdx(N-1);       % net flux out of last box 
                               % B.C. allows model to drain
                            
    % constant rainfall
    
    rain_factor = 0.000003;
              R = rain_factor.*ones(1,N); 
 
    % infiltration
    
    I(1:N) = 0*ones(1,N);    % constant infiltration along slope
 
    
    % water height changes
    
    dHdt = -dQdx + R - I;    % conservation statement 
    H    = H + dHdt.*dt;     % update total height of water
    Hneg = find(H<0);        % find negative water heights
    H(Hneg)=0;               % set negative water heights to zero
    
    % noise filter
    
    windowsize = 11;
    b = (1/windowsize)*ones(1,windowsize);
    a = 1;
    smoothH = filter(b,a,H);
    
    % fix H edges
    
    jmax = length(Hedge);
    for j=1:jmax
        Hedge(j) = (smoothH(j)+smoothH(j+1))/2; 
        
        % this calculates thickness of water at edges by averaging boxes
        
    end 
    
    % update elevations
    
        waterelevation = earthelevation + H;
    
    if (rem(t(i),tplot)==0)
       
        % water thickness plot
        
            subplot(1,3,1)
            plot(x,smoothH*100, 'k')
            yy = [bottomline, fliplr(smoothH*100)];
            fill(xx,yy,[.2 .3 .8]);
 
        % plot formatting
            
            title('flow thickness over time');
            ylabel('flow thickness (cm)');
            xlabel('distance downslope (m)');
            set(gca,'fontsize',18,'fontname','arial')
            ht=text(300,2,['  ',num2str(floor(round(t(i))/60)),...
                ' minutes  '],'fontsize',18); % print time in animation
            axis([0 xmax 0 2.5])
        
        % flux plot
        
            subplot(1,3,2)
            plot(x,Q*10000,'Color',[.2 .3 .8],'linewidth', 4)
            str = {'discharge per unit width (cm^2/s)'};         
            title('discharge over time');
            ylabel(str);
            xlabel('distance downslope (m)');
            set(gca,'fontsize',18,'fontname','arial')
            axis([0 xmax 0 65])
                     
         % hydrograph plot
         
            subplot(1,3,3)
            plot(t/3600,Hydro_1*10000,'r','linewidth', 3)
            hold on
            plot(t/3600,Hydro_2*10000,'b','linewidth', 3)
            plot(t/3600,Hydro_3*10000,'g','linewidth', 3)
            title('hydrograph monitoring')
            xlabel('time (hours)');
            legend('bottom of slope','middle of slope','top of slope',...
            'location','northwest')
            ylabel(str);
            set(gca,'fontsize',18,'fontname','arial')
            axis([0 tmax/3600 0 65])
            
        pause(0.05)
        hold off
        
    end 
    
end

