%% NEW CODE FOR CORAL GROWTH MODEL
%% COLE C. PAZAR GEOL 5700
%%  UPDATED: 16:25:20 ON 2.2.2016

figure(1)
clf
clearvars
       
step = 200;   %  essentially this is chooses your matrix size
G_m  = 0.01; %  initial growth rate (m/yr)
I_o  = 2000;  %  surface light intensity (µE m^-2 s^-1)
I_k  = 250;   %  saturating light intensity (µE m^-2 s^-1)
k    = 0.1;   %  extinction coefficent (m^-1)

% depth array (m)
Z_min = 0;
Z_max = 500;
dZ    = Z_max/step;
Z_i   = Z_min:dZ:Z_max-1;

% distance array (m)
X_min = 0;
X_sub = 100000; % 50 km to subduction zone from shoreline
X_max = X_sub;
dX    = X_max/step;
X     = X_min:dX:X_max-1;

% temporal array (ka)
years = 12000; % how many years do you want?  12 ka
t_max = years;
dt    = 10;
t     = 0:dt:t_max-1;

% plotting controls
imax   = length(t);
nplots = 15;
tplot  = t_max/nplots;
nframe = 0;

% flexure code
rho_m = 3300;
rho_w = 1030;
g     = 9.8;
D     = 10^22;
alpha = ((4*D)/((rho_m-rho_w)*g))^(1/4);
Wo    = 10000; % maximum deflection = 10 km
c     = 0.1;   % convergence rate [=] m yr^-1 

A_SL  = 60;   % 120 meters of oscillation 
P     = 22000; % 22 ka period 

beach = 10;                 %  10 meters deep water near the beach !
a     = 5*10^-8;            %  topographic coefficent [=] 1/hillslope diff.
Z_b  = (X.^2)*a + beach;   %  bed topography eqaution diffusive, parabolic
% Z_b  = (Wo).*exp((-X)/(alpha)).*cos((-X)/(alpha))+10;
% delta oxygen 18 curve for the past 2 million years
dHdt_i = G_m*(tanh(I_o*exp(-k.*Z_i)/I_k)); % steady state growth rate

% time_ka = Pleist_oxy(:,1);    % grabs the data from column 1
% delta18O = Pleist_oxy(:,2);   % grabs the data from column 2
% uncertainty = Pleist_oxy(:,3);

% Z_18 = % ocean depth based on the delta O 18 curve

%% run 

for i = 1:imax
    
    % W = (Wo).*exp((-X)/(alpha)).*cos((-X)/(alpha))+1000;                
        
    dWdt = (c/alpha)*Wo.*exp(-c*(t)/(alpha)).*(-cos(c*(t)/(alpha))...
                                          -sin(c*(t)/(alpha)));
    % Z_SL = zeros(1,length(X))+A_SL*sin(2*pi*t(i)/P); % if you want an
                                                       % oscillatory 
                                                       % sea level
      Z_SL = zeros(1,length(X));

    depth = Z_b-Z_SL; % a changing depth
                                        
    dHdt = G_m*(tanh(I_o*exp(-k.*depth)/I_k)); % steady state growth rate

    H = (dHdt*t(i)); % coral heights

    S = 0.005; % flexural subsidence rate (10 mm/year)
    % Z = (Z_b+S*t(i)) - H; % update Z's with subsidence
    Z = Z_b - H; % update Z's without subsidence
    
    subareial = find(Z<0);
    Z(subareial) = 0;
    
if rem(t(i),tplot)==0
   nframe=nframe+1;
   figure(1)
    subplot('position',[0.05 0.1 0.7 0.85])
    plot(X,Z,'linewidth',2)   % corals
    pause(0.01)
       set(gca,'YDIR','reverse','fontsize',18,'fontname','arial')
       title('coral reef growth model over 120 ka','fontname','arial',...
           'fontsize',18)
       xlabel('distance to trench (km)','fontsize',18,'fontname','arial')
       ylabel('depth (m)','fontname','arial','fontsize',18)
       axis([0 100000 -50 500])
       grid on
    hold on
    plot(X,Z_SL,'b','linewidth',1.5)      % sea level Z = 0
    plot(X,Z_b,'k','linewidth',1)         % bed topography
    legend('(colors) coral growth lines','sea level',...
        'ocean floor topography')
     subplot('position',[0.8 0.1 0.15 0.85])
      plot((dHdt*1000),Z_i,'r','linewidth',1.5) % growth rate vs. depth
      hold on
      plot(X-2,zeros(length(X)),'b','linewidth',1.5)   % sea level
        set(gca,'YDIR','reverse','fontsize',18,'fontname','arial')
        set(gca,'yaxislocation','right')
        legend('10tanh(I/I_k)e^-^k^z')
        title('growth rate vs. depth','fontname','arial','fontsize',18)
        xlabel('dH/dt (mm/yr)','fontsize',18,'fontname','arial')
        axis([-2 12 -50 500])
        grid on
end
     
end

%% finalize

