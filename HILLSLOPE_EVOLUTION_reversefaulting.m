%% CODE FOR HILLSLOPE EVOLUTION WITH A FAULTING LANDSCAPE
%% WRITTEN BY COLE C. PAZAR, ON FEBURARY 8TH, 2016.
%  1D diffusion staggered variables code style inspired by RSA

%  for a reverse fault now

%  clear all

clearvars
figure(1)
clf

%% initialize

k = 0.2;     % hillslope diffusivity [=] m^2/year

% distance array [=] m

xmin = -2000;
xmax = 2000;
dx = 1;
x = xmin:dx:xmax;

% initial topography and soil thickness [=] m
disp = 1000;
slopedegrees = -45;      % slope 
m = tand(slopedegrees);  % slope re-calculated
zmid = 1000;              % elevation at middle of an already eroded hillslope
% zb = (m.*x) + zmid;      % topography
zb = -((0.0001).*(x+disp).^2) + zmid;  % topography

H = zeros(1,length(x));  % soil thickness
z = zb+H;                % net
topo = -((0.0001).*(x+disp).^2) + zmid;
% weathering and transport parameters

wnot = 10^-3; % [=] m/yr
wstar = 0.2;  % scale for weathering rate, 1 meter
A = 0.01;
hstar = 0.2;
wdot = zeros(1,length(x)); % sets up array for weathering rate

% channel incision

edot = 10^-3; % [=] m/yr

% temporal array

tmax = 110000;
dt = 1;
t = 0:dt:tmax;
nplots = 110;
tplot = tmax/nplots;

% faulting
random = abs(randn(1,length(x))); % trying to incorporate randomness
dipdegrees = 45;
dip = tand(dipdegrees);
fault = (dip.*(x+300))+(zmid-500);
horzslip = 100/sqrt(2);           
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
    H = max(0,H);
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
    %subplot(3,1,1)
    subplot('position',[0.08 0.4 0.85 0.55])
 %  subplot('position',[left bottom width height])
    plot(x,topo,'k--','linewidth',1.5) % initial condition
    hold on
    plot(x,z+2,'g','linewidth',4)
    hold on
    plot(x,zb,'k','linewidth',2)
    hold on
    plot(x(faultzone),fault(faultzone),'r-.','linewidth',2)
    ht=text(450,900,['  ',num2str(t(i)/(1000)), ' ka '],'fontsize',18); % this makes the years print on screen
    axis([-2000 2000 0 2000]) 
    title('1D landscape evolution model of a 40º reverse-faulted, initially parabolic hillslope, 100 m of slip every 20 ka','fontsize',20)
    xlabel('distance [m]','fontsize',18)
    ylabel('elevation [m]','fontsize',18)
    set(gca,'fontsize',18)
    legend('initial condition','regolith','topography','fault plane','Location','northeast')

    hold off
    
    %subplot(3,1,3)
    subplot('position',[0.08 0.08 0.85 0.25])
 %  subplot('position',[left bottom width height])
    plot(x,H+2,'g','linewidth',4)
 %  axis([-2000 2000 0 15])
    ylabel('regolith thickness [m]','fontsize',18)
    set(gca,'fontsize',18)
    grid on
    
    drawnow
    
    end
   
end

%% end of code
