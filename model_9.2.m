%% code to simmulate the ordering of sand tracjetories into ripples
% 1d ripples vanilla

clear global
clearvars
figure(1)
clf

%% parameters used

step = 200; % choose the array lengths
font = 16;  % choses font size
d = 0.001; % grain size of the sand
porosity = 0.30; % porosity of the sand

%% initialize arrays

% x array
bin = 10*d; % bin size, also known as 'dx'
x_min = 0;
x_max = 1; % meter
x = x_min:bin:x_max;

% z array
z = zeros(size(x)); % zero and linear initial topogrpahy
bins = (z*bin)/(d^2); 

% time array
dt = 0.1; 
tmax = 1000; 
t = 0:dt:tmax;

% particle trajectory
alpha = 10; % impact angle
tan_alpha = tand(alpha); % tangent of impact angle
position = ((tan_alpha)*x_max);
N = 10; % rate of the gun's firing, basically a correction variable
splash = 4; % how much splashing happens?
v_splash = 0.1;  % splash velocity? of 1/10th of a meter per second

% plotting controls
i_max = length(x);
n_plots = step; % how many plots?
t_plot = tmax/n_plots;

%% run loops
for i = 1:length(t) 
    flux=(tan_alpha*x)+rand(1)*position;
         sub=find(flux<z);
    if isempty(sub)==1
        impact = 1+floor(i_max*rand(1));
    else
        impact = sub(1); 
    end    
    bins(impact) = bins(impact)-splash;
       bins(max(1,mod(impact+splash,i_max)))=bins(max(1,mod(impact+splash,i_max)))+splash; 
    z = N*bins*(d^2)/(4*(1-porosity)*bin); % update z
        z(1) = z(end); % fixes the B.C.
if  (rem(t(i),t_plot)==0)
        figure(1)
        fill(x*100,z*100,'k','linewidth',2); % plot topography
        axis([0 x_max*100 0 10])
        grid on
        xlabel('distance, x [cm]','fontsize',font,'fontname','arial')
        ylabel('elevation, z [cm]','fontsize',font,'fontname','arial')
        set(gca,'fontsize',font,'fontname','arial')
        title('1d ripple formation')
        drawnow
        pause(0.01)
        
end
end
