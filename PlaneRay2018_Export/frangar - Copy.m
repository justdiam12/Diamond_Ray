% frangar.m
function [c a a1 a2 a3] = frangar(T,S,D, ph,f)
%   
% Calculates sea water acoustic attenuation acording to Francois and Garrison
% J.Acoust.Soc Am. vol 72 no 6 November 1982.
% Input values
%     T = temperature  in degree centigrades 
%     S = salinity in ppm
%     D =Depth in meter
%     c= sound speed in m/s
%     ph = pH value
%     f = frequency in kHz.
%Result is
%  
%   a1= boric acid contribution to the absorption in dB/m
%   a2= Magnesium sulfate contribution to the absorption in dB/m
%   a3=pure water  contribution to the absorption in dB/m
%   a= Total absorption in sB/m

theta=T+ 273;
c=1412+3.21*T+ 1.19*S +0.0167*D;


% Boric acid contribution.
A1= (8.86/c)* 10^( 0.78*ph-5);
P1=1;
f1= 2.8*((S/35)^0.5)*10^(4-1245/theta); 
% Magnesium sulfat contribution
A2= 21.44 *(S/c)*(1+0.025*T);
P2=1-1.37*10^(-4)*D+6.2*10^(-9)*D^2;
f2=(8.17*10^(8-1990/theta))/(1+0.0018*(S-35));

%Pure water contribution
if T<=20; A3=4.937*10^(-4)-2.59*10^(-5)*T+ 9.11*10^(-7)*T^2-1.5*10^(-8)*T^3;end;
if T>20;  A3=3.964*10^(-4)-1.46*10^(-5)*T+ 1.45*10^(-7)*T^2-6.5*10^(-10)*T^3;end;

P3=1-3.83*10^(-5)*D+ 4.9*10^(-10)*D^2;

temp1= (A1.*P1.*f1.*f.^2)./(f.^2+ f1.^2); 
temp2= (A2.*P2.*f2.*f.^2)./(f.^2+ f2.^2);
temp3= A3*P3.*f.^2;
a1=temp1;a2=temp2;a3=temp3;

a = temp1+ temp2+ temp3;


