%% CompareRef
%% Jens M. Hovem august 2017;
 
 
%% Model parameters
clear ;close all;
ap=0; as=0;
c_0 =1500;c_p1=1700;rho_0=1000;rho_1=1500;rho_2=1500;
c_p2=3000;c_s1=0;c_s2=500;layer=1;
rho0=1000; ap0=.5;
a_p1=0.5;a_p2=0.5;a_s1=0; a_s2=0.5;

%a_p1=ap;a_p2=ap;a_s1=as; a_s2=as;
D=100; zs=50;z=99; eps=10^-4;
% Pekeris example 17_4_6_1
ap=0; as=0;
c_0 =1500;c_p1=1700;rho_0=1000;rho_1=1500;rho_2=rho_1;
c_p2=c_p1;c_s1=eps;c_s2=eps;layer=0;
rho0=1000; ap0=.0;
a_p1=0.5;a_p2=0.5;a_s1=0; a_s2=0.5;
%Conversion to complex speeds;
konst=8.686*2*pi;
d_p1=a_p1/konst;d_p2=a_p2/konst;d_s2=a_s2/konst;
c_p1=c_p1*(1+1i*d_p1);c_p2=c_p2*(1+1i*d_p2);c_s2=c_s2*(1+1i*d_s2);


k=linspace(0.00001,1,10000);
frequency=[100];
for n=1:length(frequency)
    
    f=frequency(n);
    omega=2*pi*f;     k_max=omega/c_0;
    theta_0=acos(k*c_0/omega);%
    theta_deg=theta_0*180/pi;
    R3(n,:)=R_bottom(real(theta_0),f,c_0,c_p1,c_p2,c_s2,rho_1,rho_2,layer);
    R2(n,:)=R_thin_layer(real(theta_0),f,c_0,c_p1,c_p2,c_s2,rho_1,rho_2,layer);
    R1(n,:)=R_layer_k(k,f,c_0,c_p1,c_p2,c_s2,rho_1,rho_2,layer);
    
end

Refloss1= -20*log10(abs(R1));
Refloss2= -20*log10(abs(R2));
Refloss3= -20*log10(abs(R3));
ang=real(theta_deg);
figure(1); %Reflection loss as function of angle

fonts= 16; lines =2;
set (gcf, 'DefaultTextFontSize', fonts, 'DefaultAxesFontSize', fonts);
plot(ang,Refloss1,'r',ang, Refloss2,'b',ang,Refloss3, 'g--','linewidth',lines);axis ij;grid;
axis([ 0 90 0 40]); xlabel( 'Grazing angle - deg');
ylabel('Reflection loss - dB')
legend (' R1',' R2',' R3')
xmin=0.2;ymin =1;dy=-3;
text(xmin, ymin-dy, ['Freq.= ' num2str(f) ' Hz']);
text(xmin, ymin-2*dy, [' D =   ' num2str(D) ' m']);
text(xmin, ymin-3*dy, [' SD =   ' num2str(zs) ' m']);
text(xmin, ymin-4*dy, [' RD =   ' num2str(z) ' m']);

figure(2);
fonts= 16; lines =2;
set (gcf, 'DefaultTextFontSize', fonts, 'DefaultAxesFontSize', fonts);
plot(ang,abs(R1),'r',ang,abs(R2),'b', ang,abs(R3),'g--', 'linewidth',lines);
axis([ 0 90 0 1.1]); xlabel( 'angle -deg');grid;
ylabel('Reflection coeff')
legend (' R1',' R2',' R3')
xmin=0.2;ymin =1;dy=0.1;
text(xmin, ymin-dy, ['Freq.= ' num2str(f) ' Hz']);
text(xmin, ymin-2*dy, [' D =   ' num2str(D) ' m']);
text(xmin, ymin-3*dy, [' SD =   ' num2str(zs) ' m']);
text(xmin, ymin-4*dy, [' RD =   ' num2str(z) ' m']);

figure(3);
fonts= 16; lines =2;
set (gcf, 'DefaultTextFontSize', fonts, 'DefaultAxesFontSize', fonts);
plot(k ,abs(R1),'linewidth',lines);
xlabel( 'Wave number 1/m');    ylabel('Reflection coeff')
legend (' R1')
xmin=0;ymin =1;dy=0.05;
text(xmin, ymin-dy, ['Freq.= ' num2str(f) ' Hz']);
text(xmin, ymin-2*dy, [' D =   ' num2str(D) ' m']);
text(xmin, ymin-3*dy, [' SD =   ' num2str(zs) ' m']);
text(xmin, ymin-4*dy, [' RD =   ' num2str(z) ' m']);

%
%
%
% Bottom reflection loss multi-frequency plot

theta0=0.01:0.01:pi/2;
fmin=0; fmax=2500;
delta_f=1;
N_freq=fix((fmax-fmin)/delta_f);



for n=1:N_freq
    f=1*n;
    F(n)=f;
    RR(n,:)=R_bottom(theta0,f,c_0,c_p1,c_p2,c_s2,rho_1,rho_2,layer);

end;

    Rabs=abs(RR);
    Rabs(Rabs<0.1)=0.1;
    RL=-20*log10(Rabs);




ma=max(max(RL)); mi=min(min(RL));
%ma=40;
mi=0; ma=10;
mst=(ma-mi)/30;
cvec=[mi:mst:ma];
dr=180/pi;
figure(54);clf; fonts=16; lines=2;
set (gcf, 'DefaultAxesFontSize', fonts );set (gcf, 'DefaultTextFontSize', 12 );
contourf(theta0*dr,F,RL, cvec)
%contourf(theta0*dr,F,RL)
colormap(flipud(gray));colorbar;grid
%title (' Reflection loss – dB' );
xlabel('{Grazing angle – deg }')
ylabel('{X = d·f - m/s }')
legendtext=0;
if legendtext==1;
    numc0=num2str(c_bottom);
    numcp1=num2str(real(cp1));
    numcp2=num2str(real(cp2));
    numcs2=num2str(real(cs2));
    numrho0=num2str(rho0);
    numrho1=num2str(rho1);
    numrho2=num2str(rho2);
    df=fmax/15;
    f1=df; f2=2*df; f3=3*df;f4=4*df;f5=5*df; f6=6*df;
    g1=2.5;g2=10;g3=20;
    text(g1,fmax-f1,['{c_0 = } ' num2str(c_bottom),' m/s'])
    text(g1,fmax-f2,['{c_{p1}=} ' num2str(real(cp1)) ' m/s'])
    text(g1,fmax-f3,['{c_{p2}=} ' num2str(real(cp2)) ' m/s'])
    text(g1,fmax-f4,['{c_{s2}=} ' num2str(real(cs2)) ' m/s'])
    text(g1,fmax-f5,['{\rho_1=} ' num2str(real(rho1)) ' kg/m^3'])
    text(g1,fmax-f6,['{\rho_2=} ' num2str(real(rho2)) ' kg/m^3'])

end


%