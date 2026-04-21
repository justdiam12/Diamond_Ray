

% PLR_initial_settings
% Default parameters that can be changed for the case you consider 

range_source=0; % Source position in range
islinear=0;     % When 0 use linear gradients, circular raypaths
frequency=[ 25 50 100 200 ];

Ne=1;% Number of elements in the source array
delta_x=0; %Element spacing in the source array
carrier_frequency=0;
lay_thick=0; % The sediment layer thickness nmust be constant, 
sigma_surface=0;% Sea surface roughess parameter in meter
sigma_bottom=0;% Bottom surface roughess parameter in meter
fs=1024;
nfft=2048;
conplot=0;
tlplot=1;
timeplot=1;

c_red=1500;
cp1=1700;
ap1=1.0;
cp2=1700;
ap2=0.1;
cs2=0.0001;
as2=0.1;
rho1=1500;
rho2=1500;
temperature=5;
salinity=35;
bothit=100;
surface_stop=0;
bottom_stop=0;
sspmethod='spline';
bathymethod='pchip';

N=2000000;
countlength=10000;

searchrays=0;
diagnostics=0;
save_in_memory=1;
eigrange=[];
eigplot=0;
envplot=1;
quiet=0;

eigprecision=[];
range_receiver=[];
range_phone=[];
z_max=[];