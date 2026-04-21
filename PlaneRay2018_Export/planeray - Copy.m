% %   Author: Jens M. Hovem
%   Date: 01:07:2016

%% planeray.m
%% The calculations start by calling the script file planeray.m
%% This program calls on the all the relevant PlaneRay programs and the input
%% programs planerayinput to set up the environment description and the
%% parameters for the calculations stored in the struct files env. and para.
close all; clear;

%% Set up the parameter files env. and para struct files containing all
%% environmental parameters and parameters for the calculations.
Example=input( 'Example ? ');
rawinput = planerayinput(Example);
[env, para] = initpara(rawinput);
save env; save para;

%% Displays the content of the env. and para. files
%% If this scenario is not what you want, stop at this point and go back.
if ~para.quiet;
    disp('Environmental parameters');disp(env);
    disp('Run parameters'); disp(para);
end
[HIST, COUNT, SUFBOT, Rays] = tracerays(env, para);

%% Main program for sorting and interpolation of the ray history
%% and produces the HIST results.
[eigenangle, HIST, COUNT, SUFBOT, Rays, para]=...........
    sortrays(HIST, COUNT, SUFBOT, Rays, env, para);

%% Calculate and plot the transmission loss for the frequencies specified
%% in para.frequency in figure 10 (and 11 for the contour plot)
if para.tlplot==1;TL=transloss(HIST, COUNT, SUFBOT, env, para, -1);end;

%% Time domain results.
para.timeplot=1;

if para.timeplot==1;
    
    %% The source pulse is generated in the getsourcesignal(source_type,para)
    %% With source_type=1, the source signal is a Ricker pulse.
    %% Other source signals can be implemented by the user
    %% The time responses at ranges specified in para.range_phone; are produced
    type=1;    [source_signal, t_start] = getsourcesignal(type,para);
    [signal, h, eigen]=timeresponse(HIST, COUNT, SUFBOT, env, para, source_signal, t_start);
    plot_time_signals;
    %% The received signals are plotted in various forms as function of both real and reduced times  ;
    %% The time plots are numbered from 20 to 24
end    
 %%% ******************* Validation and diagnostics*********************
 
if para.diagnostics==1
    threshold=0; plots=55;
    [ eigresult] = validate( para, env, SUFBOT,eigen, threshold, plots);
    % Plots only eigen rays with amplitude greater than the value of threshold
    ploteigstructure(Rays, HIST, COUNT, env, para)
    plotting=0;
    [caustics ] = FindCaustics(HIST, eigen,para, env, 200 );
    
end%

%% Plots  the bottom and surface reflection losses as function of angle
%% and frequency
plotrefloss=0;
if  plotrefloss>0;plot_reflectionloss;end

if para.timeplot==1;%% Plots the source sgnal and its frequency spectrum in figure 50
    source_sgnal=1; figure_no=50;
    plot_source_signal(source_signal,para, figure_no);
end

%%Plots source directivity  in figure 51
if para.Ne>1;
    % Specify the element spacing in mete for instance delta_x=10; 
    plot_beamspectrum;
end
    


 

