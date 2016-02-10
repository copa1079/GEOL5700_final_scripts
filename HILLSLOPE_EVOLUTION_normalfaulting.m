%% CODE FOR HILLSLOPE EVOLUTION WITH A FAULTING LANDSCAPE
%% WRITTEN BY COLE C. PAZAR, ON FEBURARY 8TH, 2016.
%  1D diffusion staggered variables code style inspired by RSA

%  normal faulting !!!

%  clear all

clearvars
figure(1)
clf

%% initialize

k = 0.1;     % hillslope diffusivity [=] m^2/year

% distance array [=] m

xmin = -2000;
xmax = 2000;
dx = 1;
x = xmin:dx:xmax;

% initial topography and soil thickness [=] m

slopedegrees = -25;      % slope 
m = tand(slopedegrees);  % slope re-calculated
zmid = 5000-40;              % elevation at middle of an already eroded hillslope
% zb = (m.*x) + zmid;      % topography
zb = -((k/100).*(x+100).^2) + zmid;  % topography
topo = -((k/100).*(x+100).^2) + zmid;  % topography

H = zeros(1,length(x));  % soil thickness
z = zb+H;                % net

% weathering and transport parameters

wnot = 10^-3; % [=] m/yr
wstar = 1;  % scale for weathering rate, 1 meter
A = 0.01;
hstar = 1;
wdot = zeros(1,length(x)); % sets up array for weathering rate

% channel incision

edot = 10^-2; % [=] m/yr

% temporal array

tmax = 188000;
dt = 1;
t = 0:dt:tmax;
nplots = tmax/1000;
tplot = tmax/nplots;

% fault w/random variables

dipdegrees = -75;
dip = tand(dipdegrees);
fault = (dip.*(x+300))+zmid-2000;
horzslip = 100;              % on average 20 meters of slip every 1000 years
vertslip = horzslip*(tand(dipdegrees)-tand(slopedegrees)); % accounts for geometry of the fault
slipperiod = 20000;          % fault period: every 20 ka


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
    
    % now set up the faulting
    faultzone = find(fault<=zb);
    if t(i)~=0 & (rem(t(i),slipperiod)==0) % remove the (i) on slipperiod to fix?
     % faultstart=find(fault<=zb);
     % if t(i)==20000
        slip=ones(size(faultzone)).*vertslip;
    
        for j=1:length(faultzone)
            if (zb(faultzone(j))+slip(j))<=fault(faultzone(j))
                zb(faultzone(j))=fault(faultzone(j));
                z(faultzone(j))=fault(faultzone(j));
                H(faultzone(j))=0;
            
            else
                zb(faultzone(j))=zb(faultzone(j))+slip(j);
                z(faultzone(j))=z(faultzone(j))+slip(j);  
                
            end
            
        end
      
        end
   
    if(rem(t(i),tplot)==0)
        
    figure(1)
    %subplot(2,1,1)
    subplot('position',[0.08 0.4 0.85 0.55])
    plot(x,topo,'r--','linewidth',1.5)
    hold on
    plot(x,z+40,'g','linewidth',3)
    hold on
    plot(x,zb,'k','linewidth',2)
    hold on
    plot(x(faultzone),fault(faultzone),'c-.','linewidth',2)
    ht=text(-1800,4250,['  ',num2str(t(i)/1000), ' ka '],'fontsize',18); % this makes the years print on screen
    axis([-2000 2000 -400 5400]) 
    title('Landscape evoltuion of a 75º normal-faulted hillslope, 100 m of horizonatal slip every 20 ka','fontsize',16)
    xlabel('Distance [m]','fontsize',18)
    ylabel('Elevation [m]','fontsize',18)
    set(gca,'fontsize',16)
    
    legend('Initial profile','Regolith','Bedrock','Fault plane','Location','northeast')

    hold off
    
    %subplot(2,1,2)
    subplot('position',[0.08 0.08 0.85 0.25])
    plot(x,H,'g','linewidth',3)
    % axis([-2000 2000 0 30])
    % title('Regolith depth as a function of distance','fontsize',16)
    ylabel('Regolith thickness [m]','fontsize',18)
    xlabel('Distance [m]','fontsize',18)
    set(gca,'fontsize',16)
    grid on
    
%     subplot(3,1,3)
%     plot(x,H*100,'g','linewidth',3)
%     axis([-2500 2500 0 100])
%     title('regolith depth as a function of distance','fontsize',20)
%     xlabel('distance [m]','fontsize',15)
%     ylabel('regolith depth [cm]','fontsize',15)
%     set(gca,'fontsize',15)
%     grid on
    
    drawnow
    
    end
   
end

%% end of code
