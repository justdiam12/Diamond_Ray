function [env, para] = initpara(rawinput)
%initpara

%   Authors: Jens M. Hovem and Shefeng Yan
%   $Revised 7. April 2014 $ filtering of ang_slope;

disp(rawinput.title_tekst);

if isempty(rawinput.z_max),
    rawinput.z_max=max(rawinput.z_input);
end;
if ~strcmpi(rawinput.sspmethod,'none')
    z_interpolated = 0:rawinput.del_z:rawinput.z_max;
    c_interpolated = interp1(rawinput.z_input.',rawinput.c_input.',z_interpolated,rawinput.sspmethod);
else
    z_interpolated = rawinput.z_input;
    c_interpolated = rawinput.c_input;
end;
if ~all(diff(c_interpolated))
    rawinput.rayopt=1;
else
    rawinput.rayopt=0;
end;
ssp(:,1)           = [z_interpolated.';z_interpolated(end)+z_interpolated(2)-z_interpolated(1)];
ssp(:,2)           = [c_interpolated.';c_interpolated(end)];
[R_bathy, Z_bathy, R_slope, ang_slope] = initbathymetry(rawinput.bathymethod,rawinput.R_b, rawinput.Z_b, rawinput.del_z);
N_filter=16;
A=1; B=ones(1, N_filter)/N_filter;
ang_slope=filter(B, A, ang_slope);
if isempty(rawinput.range_receiver)
    range_receiver     = linspace(0,rawinput.Rmax,401);
    rawinput.range_receiver=range_receiver(2:end);
end;
if isempty(rawinput.range_phone)
    rawinput.range_phone=rawinput.range_receiver;
end;


start_theta= sort(rawinput.start_theta);
start_theta= start_theta(start_theta>-90 & start_theta<90);
start_theta=start_theta(start_theta~=0);
start_theta(abs(diff(start_theta))<=10*eps)=[];
rawinput.start_theta=start_theta;
if length(rawinput.start_theta)>301
    rawinput.savememory=1;
    disp('Number of angles larger than 301, and the ray trace data are not saved.');
    disp('savememory = 1.');
end;
if length(rawinput.laythick)==1 || ~isfield(rawinput,'RR')
    rawinput.RR=[0,rawinput.Rmax]; 
end;
if length(rawinput.laythick)>=length(rawinput.RR)
    rawinput.laythick=rawinput.laythick(1:length(rawinput.RR)-1);
end;

para.del_z         = rawinput.del_z;
para.N             = rawinput.N;
para.Rmax          = rawinput.Rmax;
para.start_angles          = rawinput.start_theta;
para.range_source  = rawinput.range_source;
para.z_source      = ceil(rawinput.z_source/rawinput.del_z)*rawinput.del_z;
para.z_receiver    = ceil(rawinput.z_receiver/rawinput.del_z)*rawinput.del_z;
para.range_receiver= rawinput.range_receiver;
para.range_phone   = rawinput.range_phone;
para.eigrange      = rawinput.eigrange;
para.frequency     = rawinput.frequency;
para.fs            = rawinput.fs;
para.nfft          = rawinput.nfft;
para.conplot       = rawinput.conplot;
para.tlplot        = rawinput.tlplot;
para.timeplot      = rawinput.timeplot; 
para.envplot       = rawinput.envplot;
para.eigplot       = rawinput.eigplot;
para.rayopt        = rawinput.rayopt;
para.c_red         = rawinput.c_red;
para.bottom_stop   = rawinput.bottom_stop;
para.surface_stop  = rawinput.surface_stop;
para.Ne            = rawinput.Ne;
para.delta_x       = rawinput.delta_x;
para.bothit        = rawinput.bothit;
para.carrier_frequency= rawinput.carrier_frequency;
para.countlength   = rawinput.countlength;
para.diagnostics   = rawinput.diagnostics;
para.quiet         = rawinput.quiet;
para.bathymethod = rawinput.bathymethod;
para.sspmethod = rawinput.sspmethod;

env.title          = rawinput.title_tekst;
env.ssp            = ssp;
env.bathy          = [R_bathy; Z_bathy];
env.bathslope      = [R_slope; ang_slope];
env.sigma_surface = rawinput.sigma_surface; %Sea surface roughness
env.sigma_bottom = rawinput.sigma_bottom; %Bottom surface roughness
env.RR       = rawinput.RR;
env.cp1      = rawinput.cp1;
env.ap1      = rawinput.ap1;
env.rho1     = rawinput.rho1;
env.cp2      = rawinput.cp2;
env.cs2      = rawinput.cs2;
env.ap2      = rawinput.ap2;
env.as2      = rawinput.as2;
env.rho2     = rawinput.rho2;
env.laythick = rawinput.laythick;
env.temperature=rawinput.temperature;
env.salinity=rawinput.salinity;



%% initbathymetry
function [R_bathy, Z_bathy, R_slope, ang_slope] = initbathymetry( bathymethod, R_b, Z_b, del_z)
% Modification 4 April 2014
Nb=401;
del_z=del_z/2;
newR_b=linspace(min(R_b),max(R_b),Nb); 
newZ_b=interp1(R_b,Z_b,newR_b, bathymethod, 'extrap');
R_b=newR_b;
Z_b=newZ_b;
Z_b=round(Z_b/del_z)*del_z;
kk=round(Z_b/del_z);

sameind1=find(diff(kk)==0)+1;
sameind2=find(diff(sameind1)==1);
R_b(sameind1(sameind2))=[];
Z_b(sameind1(sameind2))=[];
kk=round(Z_b/del_z);

Z_bathy=[];
R_bathy=[];
omit=0;
for ii=1:length(kk)-1
    if kk(ii)<kk(ii+1)
        if ii==length(kk)-1
            thisz  =(kk(ii):kk(ii+1))*del_z;
        else
            thisz  =(kk(ii):kk(ii+1)-1)*del_z;
        end;
        thisr  =interp1([Z_b(ii),Z_b(ii+1)], [R_b(ii),R_b(ii+1)], thisz, bathymethod, 'extrap');
        if omit
            thisz=thisz(1+omit:end);
            thisr=thisr(1+omit:end);
            omit=0;
        end;
    elseif kk(ii)>kk(ii+1)
        if ii==length(kk)-1
            thisz  =(kk(ii):-1:kk(ii+1))*del_z;
        else
            thisz  =(kk(ii):-1:kk(ii+1)+1)*del_z;
        end;
        thisr  =interp1([Z_b(ii),Z_b(ii+1)], [R_b(ii),R_b(ii+1)], thisz, 'linear', 'extrap');
        if omit
            thisz=thisz(1+omit:end);
            thisr=thisr(1+omit:end);
            omit=0;
        end;
    else
        thisz  =[kk(ii),kk(ii+1)]*del_z;omit=1;
        thisr  =[R_b(ii),R_b(ii+1)];
    end;
    Z_bathy=[Z_bathy,thisz];
    R_bathy=[R_bathy,thisr];
end;

newR_b=linspace(min(R_b),max(R_b),Nb); 
newZ_b=interp1(R_b,Z_b,newR_b, bathymethod);
R_slope=newR_b;
AA=-atan(diff(newZ_b)./diff(newR_b));
ang_slope=[AA AA(end)];

