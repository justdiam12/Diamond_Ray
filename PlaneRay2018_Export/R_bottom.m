
function R=R_bottom(theta0,f,cp0,cp1,cp2,cs2,rho1,rho2,d)
% %R_layer calculates the reflection coefficient for  a sediment layer over
% solid bedrock.
% theta0= grazing angle of incoming wave, radians
% d=layer thickness in meter
% f=frequency in Hz
%cp0=sound velocity in the water
%cp0=sound velocity in the water, assumed to be 1500 m/s
%cp1=complex sound speed the sediment
%cp2=complex comp. speed in the bedrock
%cs2=complex shear speed in the bedrock
%rho0=density of water, assumed to be 1000kg/m^3
%rho1=density of the sediment
%rho2=density of the bedrock

rho0=1000.0;
c0=cp0;
omega=2*pi*f;
k=omega*cos(theta0)/c0;
gamma0 =sqrt((omega./c0 ).^2-k.*k);
gammap1=sqrt((omega./cp1).^2-k.*k);

gammap2=sqrt((omega./cp2).^2-k.*k);
gammas2=sqrt((omega./(cs2+eps)).^2-k.*k);

if imag(gammap1)>0, gammap1=conj(gammap1);end;
if imag(gammap2)>0, gammap2=conj(gammap2);end;
if imag(gammas2)>0, gammas2=conj(gammas2);end;

costheta2s=k*cs2./(omega+eps);
sintheta2s=cs2*gammas2./(omega+eps);% test1=costheta2s.^2+sintheta2s.^2;
cos2s=2*costheta2s.^2-1;
sin2s=2*sintheta2s.*costheta2s;
C=cos2s.^2;
S=sin2s.^2;
Z0 =omega.*rho0./(gamma0 +eps);
Zp1=omega.*rho1./(gammap1+eps);
Zp2=omega.*rho2./(gammap2+eps);
Zs2=omega.*rho2./(gammas2+eps);

r01=(Zp1-Z0)./(Zp1+Z0+eps);
r12=(Zp2.*C+Zs2.*S-Zp1)./(Zp2.*C+Zs2.*S+Zp1+eps);
E=exp(-2*i*gammap1*d);
R= (r01+r12.*E)./(1+r01.*r12.*E);

