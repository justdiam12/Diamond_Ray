%getbeamspec.m
%%   Author: Jens M. Hovem
%%   Copyright 2011 Acoustic Research Center, NTNU
%%   Revised: 2011/01/11

%%  Introduces source directivity in the calulation of the acoustic field
%%  arrayno =1;gives the beam pattern of a line array with Ne elements at
%%  a spacing delta x. These parameters are loaded from file para.
%%  The beamspectrum values are used in transfunc 
%%  Other beampattern can be inserted by the user
%%  The initial angle is luanch angle of the source as specified in the
%%  input file 

function [beamspectrum ] = getbeamspec(para, initial_angle, arrayno)
Nf=para.nfft;
fs=para.fs;   
Ne=para.Ne;
if arrayno==1
    delta_x=para.delta_x;
    c=1500;
    f=(0:Nf-1)*fs/Nf;
    theta=initial_angle;
    d=delta_x;
    w=pi.*f/c.*d.*cos(theta)+eps;
    F=exp(-1i.*w*(Ne-1));
    B=F.*sin(Ne*w)./(Ne*sin(w));
    B_spectrum=B(1:Nf/2);
    beamspectrum=B_spectrum;
end;


    

    
 

      
      
      
    
       