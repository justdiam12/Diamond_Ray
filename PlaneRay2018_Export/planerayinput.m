%%  planerayinput.m
%%  Author: Jens M. Hovem
%%  Copyright 2011 Acoustic Research Center, NTNU
%%  $Revised $  $Date: 2011/01/11
%% This program sets the parameters most likely to be used, but the settings can
%% be overruled by statements below for the case you want to consider;

%% Note that the parameters  "searchrays" and "beamdisp" are not acivated in
%% this version of PlaneRay
 

%% Test cases

%% case 0 Test case
%% case Pekeris 200 depth 
%% case 2 Summer
%% case 3
%% case 4
%% case 5
%% case 
function rawinput = planerayinput(Example)


PLR_initial_settings;


switch Example;
    
    case 0
        title_tekst='Test case'; %Set the name of your case
        load ssp_case_0;  ssp=ssp_case_0;% Generate or load the sound speed profile
        c_input=ssp(:,1); z_input=ssp(:,2); del_z=1;
        c_red=1490;% reduction speed for the time scale of the time plots
        
        %Source and receives positions
        R_max=5000;    % Maximum range in meter for ray calculations
        z_source=25;    % Source depth in meter
        z_receiver=100; % Receiver depths in meter
        
        % Generate or load the bathymetry
        R_b=0:1: R_max;   z_max=max(z_input);
        Z_b=0.95*z_max-R_b*0.015+ 10*sin(pi*(5*R_b/R_max));
        bathymethod='pchip';% Chose interpolation method
        
        %Geoacoustic input
        cp1=1700; ap1=0.5; rho1=1700;
        cp2=cp1; ap2=ap1; rho2=rho1;
        cs2=0; as2=.0;lay_thick=5; %Sediment layer thickness
         fs=512; nfft=4096 ; %Sampling frequency, FFT block lengt
        sigma_surface=0;%Sea surface roughness in meters
        sigma_bottom=0;%Bottom surface roughness in meters
        
        % Specifications of initial launch angles
        N_angle=50 ;  theta_max=30; n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;ang=ang(ang>0);
        
        start_theta=sort([-ang ang]);
        % Limit the numbers of bottom/surface interactions
        bottom_stop=4;% Stop after so many bottom interactions “0) don’t stop
        surface_stop=4; % Stop after so surface bottom interactions “0) don’t stop
        frequency=[  100  1000];% For  transmission loss calculations
        % carrier_frequency=19000;% Used to calculate the seawater absorption
        
        range_receiver=10:10:R_max; % For transmission loss calculations,
        range_phone=500:500:R_max; % for time domain calculations
        %plotrefloss=1;
         diagnostics=1;
        eigplot=1;eigrange=[ 4000]; % For evaluation of eigenrays
         %For more info on the calculated wave fields
       
        
        
        
        
    case 1 %
        title_tekst='Pekeris-200 Depth'
        wd=200;
        del_z=.1; z_input=0:del_z:wd;
        c_input=1500*ones(size(z_input));N=100000;
        R_max=5000;
        
        z_source=25; z_receiver=150;
        R_b=[0:1: R_max];Z_b=ones(1,length(R_b))*wd;
        bottom_stop=0;surface_stop=0;
        sigma_surface=0;%Sea surface roughness in meters
        sigma_bottom=0;%Bottom surface roughness in meters
        N_angle=20;
        theta_max=45;n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;
        start_theta=sort([ang -ang]);
        islinear =1 ;
        fs=1024; nfft=1024;
        eigplot=1;
        
        type=input('Fluid (1)  solid (2) layered (3) ?');
        
        if type==1;title_tekst=' Fluid bottom' ;
            cp1=1700; ap1=0.5; rho1=1500;
            cp2=1700; ap2=0.5; rho2=1500;cs2=eps;as2=eps;lay_thick=0;
        end
        
        
        if type==2;title_tekst=' Elastic bottom' ;
            cp1=1700; ap1=0.2; rho1=1500;
            cp2=2500; ap2=0.2; rho2=2000;cs2=500;as2=0.50;lay_thick=0;
        end
        
        if type==3;title_tekst=' Layered bottom' ;
            cp1=1700; ap1=0.2; rho1=1500;
            cp2=2500; ap2=0.2; rho2=2000;cs2=500;as2=0.2;lay_thick=5;
            z_source=5; z_receiver=96;wd =100;
        end
               
        frequency= [ 50  100 ];
        range_receiver = 100:100:R_max;
        range_phone=0:1000:R_max;
        timeplot=1 ;
        diagnostics=1;
        
        
        
        
        
         case 2
        title_tekst='Elba E';
        load sound_speed.dat; z_input = sound_speed(:,1);c_input = sound_speed(:,2);
        load bottom.dat; rbty = bottom(:,1);zbty = bottom(:,2);
        R_b=rbty; Z_b=zbty;
        R_max=max(rbty);
        r=1:1:R_max; R_b=r;
        Z_b=interp1(rbty,zbty, r);
        z_max=max(Z_b) ;del_z=0.1;
        lay_thick=0;
        c_red=1510;
         diagnostics=1;
        z_source=120 ;z_receiver=65; range_source=0;
        bottom_stop=0;surface_stop=0;
        
        
        N_angle=100;n=1:N_angle;
        theta_max=45;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;
        ang=ang(ang>0);
        start_theta=sort([-ang ang]);
        start_theta=start_theta(start_theta~=0); % Avoid  starting a horizontal ray
        fs=1024;
        nfft=1024;
        cp1=1537; ap1=0.2; rho1=1800;
        cp2=cp1; ap2=ap1; rho2=rho1;
        cs2=0.01;as2=.0;
        frequency=[  100 ];
        carrier_frequency=25500;
        range_receiver=  0:1:R_max;
        range_phone=500:100:6500;
        sigma_bottom=0.0;
        
    case 3 %Propagation over a sea mountain
        Month=input('Summer (1) or Winter(2)');
        wd=300;
        del_z_raw=20;
        z_raw=0:del_z_raw:wd;
        del_z=.1;
        z_input=0:del_z:300;
        if Month==1;  c_raw=[1480.0 1476.7 1475.6 1475.2 1475.5 1476.5 1477.5 1478.5...
                1479.0 1479.2 1479.0 1478.8 1478.2 1478.0 1478.3 1478.6 ];
            title_tekst=('Sea Mount - Summer profile');
        end;
        
        if Month==2
            c_raw=[1472.0  1475.0 1477.0 1478.0 1478.5 1478.9 1479.0 1479.1 ...
                1479.2 1479.0 1478.8 1478.2 1478.0 1478.3 1478.6 1478.8];
            title_tekst=('Sea Mount - Winter profile');
        end;
        c_input=interp1(z_raw,c_raw,z_input,'spline');
        c=c_input;
        R_max=10000; % Maximum range in meter for ray calculations
        range_source=0;
        z_source=25; % Source depth in meter
        z_receiver=20; % Receiver depths in meter
        theta_start=-15;% Minimum angle for ray  traceing - degree
        theta_stop=5; %Maximum angle for ray  traceing - degree
        delta_theta=.1;%Increment for raytracing - angle
        start_theta=theta_start:delta_theta:theta_stop;
        
        %Definening  the botttom topograpy
        R_b=0:1000:10000;
        Z_b=[300 285 275 250 200 150 100 150 200 250 300];
        lay_thick=0;
        timeplot=0;
        range_receiver = 1:1:R_max;
        range_phone=0:1000:R_max;
        bottom_stop=1;
        surface_stop=0;
        eigrange=[  8500];
        
        
            
    case 4 %
        title_tekst='Shelf break - downslope ';
        wd=300; R_max=25000;   z_source=10; z_receiver=75;
        R_raw=[0:5000: R_max];
        r1=5000; r2=15000;
        d1=100; d2=250;
        g=(d2-d1)/(r2-r1);
        
        Z_raw=d1+g*(R_raw-r1);
        Z_raw(R_raw<r1)=d1;
        Z_raw(R_raw>r2)=d2;
        R_b= [0:1: R_max];
        Z_b=interp1(R_raw,Z_raw,R_b,'pchip');
        Z_b=fliplr(Z_b);title_tekst='Shelf break - upslope ';
        
        del_z=50; z_raw=0:del_z:wd; z1=50;
        N_z=length(z_raw);
        n=1:N_z;n0=fix(N_z*z1/wd)+1;
        g1=-0.15; g2=0.01;
        n1=1:n0; n2=n0+1:N_z;
        c1=1485+g1*z_raw(n1);
        c2=min(c1)+g2*(z_raw(n2)-z1);
        c_raw=[c1 c2];
        del_z=.1; z_input=0:del_z:wd;
        c_input=interp1(z_raw,c_raw,z_input,'spline');
        bottom_stop=12;surface_stop=20;
        start_theta=-20:0.5:20;
        fs=512;   nfft=2048;
        frequency=100; carrier_frequency=9900;
        range_phone=100:1000:R_max;
        range_receiver = 100:100:R_max;
        timeplot=1;
        cp1=1800; ap1=0.5; rho1=1500;
        cp2=2500; ap2=0.5; rho2=2500;
        cs2=500; as2=0.5;
        lay_thick=10;
        bathymethod='pchip';
        
        case 5
        title_tekst='Summer '; %Set the name of your case
        load ssp_case_0;  ssp=ssp_case_0;% Generate or load the sound speed profile
        c_input=ssp(:,1); z_input=ssp(:,2); del_z=.1;
        c_red=1490;% reduction speed for the time scale of the time plots
        
        %Source and receives positions
        R_max=3000;    % Maximum range in meter for ray calculations
        z_source=5;    % Source depth in meter
        z_receiver=100; % Receiver depths in meter
        
        % Generate or load the bathymetry
        R_b=0:1: R_max;   z_max=max(z_input);
        Z_b=0.95*z_max-R_b*0.015+ 10*sin(pi*(5*R_b/R_max));
        bathymethod='pchip';
        
        %Geoacoustic input
        cp1=1700; ap1=0.5; rho1=1700;
        cp2=3000; ap2=0.5; rho2=2500;
        cs2=0; as2=.0;
        lay_thick=19; %Sediment layer thickness
        
        sigma_surface=0;%Sea surface roughness in meters
        sigma_bottom=0;%Bottom surface roughness in meters
        
        % Specifications of initial launch angles
        N_angle=60 ;  theta_max=60; n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;ang=ang(ang>0);
        
        start_theta=sort([-ang ang]);
        ang=ang(ang>0);
        % Limit the numbers of bottom/surface interactions
        bottom_stop=2;% Stop after so many bottom interactions “0) don’t stop
        surface_stop=2; % Stop after so surface bottom interactions “0) don’t stop
        
        fs =1024; nfft=1024*4; %Sampling frequency, FFT block length
        frequency=[  50 ];% For  transmission loss calculations
        
        carrier_frequency=0;% Used to calculate the seawater absorption
        
        range_receiver =10:10:R_max; % For transmission loss calculations,
        range_phone=25:25:3000; % for time domain calculations
        
        eigplot=0;eigrange=[ 2000]; % For evaluation of eigenrays
        diagnostics=1;  %For more info on the calculated wave fields
        
          
        
    case 6
        title_tekst='Shelf break - upslope ';
        Type=input(' type = 1 for upslope and 2 for downslope');
        wd=300; R_max=50000;   z_source=6; z_receiver=15;
        R_b=[0:50: R_max];
        r1=10000; r2=40000;
        d1=100; d2=250;
        g=(d2-d1)/(r2-r1);
        Z_b=d1+g*(R_b-r1);
        Z_b(R_b<r1)=d1;
        Z_b(R_b>r2)=d2;
        
        cp1=1700; ap1=0.5; rho1=1500;
        cp2=3000; ap2=0.5; rho2=2500;cs2=600;as2=0.5;lay_thick=2;
        if Type==2;title_tekst='Upslope - elastic'; Z_b=fliplr(Z_b); lay_thick=0;end;
        if Type==1;title_tekst='Upslope - fluid'; as2=eps; cs2=eps;cp2=cp1;rho2=rho1;lay_thick=0;Z_b=fliplr(Z_b);end;
        
        del_z=1;
        z_input=0:del_z:wd;
        c_input= 1500*ones(size(z_input));
        c_red=min(c_input);
        
        
        quiet=0;
        start_theta_1=-30:5:30;
        start_theta_2=-80:5:-30;
        % start_theta=-30:5:30;
        bottom_stop=18;surface_stop=18;
        start_theta=[  start_theta_1 -start_theta_2  start_theta_2];
        islinear =0 ;
        fs=512;    nfft=2048*4;
        
        frequency=[ 25  50 100 200 ];
        carrier_frequency=0;
        range_receiver=  50:50:R_max;
        range_phone=1000:10000:R_max;
        diagnostics=1;
        timeplot=0 ; countlength=400;
        
        
     
        
    case 7 %Arctic propagation in deep water ST-586841
        load ssp_586841 ;  ssp= ssp_586841 ;
        title_tekst='ST- N 73:44.9 N;  2:40 W';
        c_raw=ssp.c; z_raw=ssp.z;
        z_trunc=z_raw(z_raw<400);   c_trunc=c_raw(z_raw<400);
        z_raw=z_trunc; c_raw=c_trunc;
        c_input=c_raw; z_input=z_raw;
        del_z=1;
        wd=300;
        R_max=10000;   z_source=295; z_receiver=100;
        R_b=[0:100: R_max];
        Z_b=ones(1,length(R_b))*wd;
        bottom_stop=1;surface_stop=0;
        
        N_angle=62; theta_max=32;n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;
        start_theta=sort([ -ang]);
        
        start_theta=start_theta(start_theta~=0);
        fs=512; nfft=1024;
        cp1=1700; ap1=0.1; rho1=1500;
        cp2=1700; ap2=0.1; rho2=1500;
        cs1=eps; as1=eps;
        cs2=0.0001; as2=0.0;lay_thick=0;
        carrier_frequency=10*10^3;
        
        temperature=0;
        salinity=25;
        
        frequency=  [ 40 60 80 110  ]*10^3
        range_receiver = 0:500:R_max;
        range_phone=0:500:R_max;
        timeplot=0;
        
        
        
    case 8 %Arctic propagation in deep water ST-604321
        load ssp604321;
        title_tekst= 'ST- 74;3N 6;25W';
        z_raw=ssp604321(:,1);       c_raw=ssp604321(:,2);
        z_trunc=z_raw(z_raw<400);   c_trunc=c_raw(z_raw<400);
        z_raw=z_trunc; c_raw=c_trunc;
        
        c_input=c_raw; z_input=z_raw;
        del_z=1;
        wd=max( z_raw);
        %wd=300;
        R_max=10000;   z_source=295 ;z_receiver=100;
        R_b=[0:100: R_max];
        Z_b=ones(1,length(R_b))*wd;
        bottom_stop=1;surface_stop=1;
        
        N_angle=30; theta_max=30;n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;
        start_theta=sort([ang -ang]);
        start_theta=start_theta(start_theta~=0);
        fs=512; nfft=4096*2;
        cp1=1700; ap1=0.1; rho1=1500;
        cp2=1700; ap2=0.1; rho2=1500;
        cs2=0.0001; as2=0.0;lay_thick=0;
        frequency=  [ 50 100 200 400  ];
        range_receiver = 0:500:R_max;
        range_phone=0:1000:R_max;
        timeplot=1 ;
        
        
    case 9
        title_tekst='73°N 44.9 2°W 40'
        load ssp_586841
        c_raw=c; z_raw=z;
        z_trunc=z_raw(z_raw<400);   c_trunc=c_raw(z_raw<400);
        z_raw=z_trunc; c_raw=c_trunc;
        c_input=c_raw; z_input=z_raw;
        del_z=1;
        wd=300;
        R_max=10000;   z_source=295; z_receiver=50;
        R_b=[0:100: R_max];
        Z_b=ones(1,length(R_b))*wd;
        bottom_stop=1;surface_stop=1;
        
        N_angle=30; theta_max=30;n=1:N_angle;
        ang=theta_max.^(n/N_angle)-theta_max.^(1/N_angle) ;
        start_theta=sort([ang -ang]);
        start_theta=start_theta(start_theta~=0);
        fs=512; nfft=4096*2;
        fs=512; nfft=1024;
        cp1=1700; ap1=0.1; rho1=1500;
        cp2=1700; ap2=0.1; rho2=1500;
        cs1=eps; as1=eps;
        cs2=0.0001; as2=0.0;lay_thick=0;
        carrier_frequency=10*10^3;
        
        temperature=0;
        salinity=25;
        
        frequency=  [ 40 60 80 110  ]*10^3
        range_receiver = 0:500:R_max;
        range_phone=0:500:R_max;
        timeplot=0;

end % case selection       
rawinput.quiet=quiet;
rawinput.countlength=countlength;
rawinput.title_tekst=title_tekst;
rawinput.del_z = del_z;
rawinput.z_input= z_input;
rawinput.c_input = c_input;
rawinput.Rmax=R_max;
rawinput.N=N;
rawinput.z_source=z_source;
rawinput.z_receiver=z_receiver;
rawinput.delta_x=delta_x;
rawinput.Ne=Ne;
rawinput.z_max=z_max;
rawinput.range_source=range_source;
rawinput.rayopt=islinear;
rawinput.R_b=R_b;
rawinput.Z_b=Z_b;
rawinput.laythick=lay_thick;
start_theta=start_theta(start_theta~=0);
rawinput.start_theta=start_theta;
rawinput.bottom_stop=bottom_stop;
rawinput.surface_stop=surface_stop;
rawinput.frequency=frequency;
rawinput.carrier_frequency=carrier_frequency;
rawinput.fs=fs;
rawinput.nfft=nfft;
rawinput.tlplot=tlplot;
rawinput.envplot=envplot;
rawinput.conplot=conplot;
rawinput.timeplot=timeplot;
rawinput.eigplot=eigplot;
rawinput.range_receiver=range_receiver;
rawinput.range_phone=range_phone;
rawinput.diagnostics=diagnostics;
rawinput.c_red=c_red;
rawinput.sigma_surface=sigma_surface;
rawinput.sigma_bottom=sigma_bottom;
rawinput.cp1=cp1;
rawinput.ap1=ap1;
rawinput.rho1=rho1;
rawinput.cp2=cp2;
rawinput.cs2=cs2;
rawinput.ap2=ap2;
rawinput.as2=as2;
rawinput.rho2=rho2;
rawinput.salinity=salinity;
rawinput.temperature=temperature;
rawinput.bothit=bothit;
rawinput.sspmethod=sspmethod;
rawinput.bathymethod=bathymethod;

rawinput.eigrange=eigrange;

rawinput.searchrays=searchrays;
rawinput.save_in_memory=save_in_memory;
rawinput.carrier_frequency=carrier_frequency;