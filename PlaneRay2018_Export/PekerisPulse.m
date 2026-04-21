%% PekerisPulse
%% Calculation of pulse response of a Pekeris waveguide
%% medium using wavenumber integration technique.
%% Jens M. Hovem July 2017;
 
 
%% Model parameters
clear ;close all;
 
 c_0 =1500;c_p1=1700;rho_0=1000;rho_1=1500;
 D=200; zs=50;z=150; eps=10^-4;
%% Setting up the Ricker pulse
f_c= 100;Nt= 1024*2;nt= 1:Nt;f_s= Nt;
delta_t= 1/f_s;time = nt*delta_t;
freq = (1:Nt)*f_s/Nt;t0 = 5e-2;
lines=2.0;fonts=18;
%% Getting the spectrum for the Ricker pulse:
S_w = zeros(1,Nt);
pulse = Ricker(f_c,t0,f_s,Nt);
spectrum = ifft(pulse,Nt);
S_w(1)= spectrum(1);S_w(2:Nt/2+1)= 2*spectrum(2:Nt/2+1);
S(Nt/2+2:Nt) = 0;
% Setting up range and horizontal wavenumber vectors:
Nr = 2048;k  = linspace(0.001,1,Nr);k_sq = k.^2;d_k  = k(2)-k(1);
d_r = 2*pi/(Nr*d_k);r = 1:d_r:2000;
% Setting up vectors and matrixes to be used in loop:
PHI=zeros(Nt,Nr);gamma_0=zeros(Nr,1);
gamma_p_1 = gamma_0;gamma_p_2 = gamma_0;gamma_s_2 = gamma_0;
%z_p_1 = zeros(Nr,1);z_p_2 = z_p_1;z_s_2 = z_p_1;
%sin_ang = zeros(Nr,1);cos_ang = sin_ang;
%r_12 = zeros(Nr,1);r_01 = r_12;

%% Calculating the wavenumber dependent field.
PHI0=zeros(Nt, Nr); PHI=zeros(Nt, Nr);  Rb=zeros(Nt, Nr);

% Perturbation of the k axis
k_real=k;k=k-i*eps;

for n=1:length(freq)
    omega=2*pi*freq(n);
    % Wavenumbers for the different media
    kappa_0 = omega/c_0; kappa_1     = omega/c_p1;
    % Vertical wavenumber for the different media
    
    gamma_0=(kappa_0^2 - k.^2).^(1/2);
    gamma_p_1 = (kappa_1^2 - k.^2).^(1/2);
 
   %Impedances for the different media
    z_p_1= omega*rho_1./gamma_p_1;
    r_01 = (gamma_0.*rho_1-gamma_p_1.*rho_0)./(gamma_0.*rho_1+gamma_p_1.*rho_0);  
    Rb(n,:)= r_01;
    % The total field: (Eq. 17.65 with Rs=-1)
    PHI0(n,:) = (1./(j*4.*pi.*gamma_0)).*(exp(i.*gamma_0.*abs(z-zs))+ Rb(n,:).*exp(i.*gamma_0.*(2.*D-zs-z)));
    PHI(n,:) = S_w(n).*PHI0(n,:);
end





%% Wavenumber integration and FFT
integrand=zeros(Nt,Nr);

figure(1);set (gcf, 'DefaultTextFontSize', fonts, 'DefaultAxesFontSize', fonts);
hold on
range=100:100:1500;
N_range=length(range);
for q=1:N_range
    select_range=range(q);
    [ a, p]=min(abs(r- select_range));
    % The range dependent part of [Eq. 17.91]
    r_dependency = (d_k/sqrt(2*pi*r(p)))*exp(i*r(p)*min(k))*exp(-i*pi/4);
    % Summation over horisontal wave number [Eq. 17.91]
    for q=1:Nr
        integrand(:,q)  = PHI(:,q)*exp(i*q*d_k*min(r))*sqrt(k(q))*exp(2*pi*i*p*q/Nr);
    end
    sum_over_k = integrand*ones(Nr,1);
    PHI_R_W = sum_over_k.*r_dependency;
    % FFT to find time response:
    PHI_R_T(:,p) = 2*real(fft(PHI_R_W))*exp(eps*r(p))*Nr;
 
    
   % Plotting the scaled pulse responses
    scale=100;
    plot(time,scale*PHI_R_T (:, p) + r(p),'b','linewidth',lines);
    grid
end
x_min=0; x_max=max(time); y_min=-100; y_max=max(range)*1.1; dy= y_max/20;
axis([ x_min x_max y_min y_max])
text(x_min, y_max-dy, ['Ricker pulse; ' num2str(f_c) ' Hz']);
text(x_min, y_max-3*dy, [' D =   ' num2str(D) ' m']);
text(x_min, y_max-4*dy, [' SD =   ' num2str(zs) ' m']);
text(x_min, y_max-5*dy, [' RD =   ' num2str(z) ' m']);

title('Pulse responses as function of range')
ylabel('Range – m');xlabel('Time – s')
box


% Plot integrand
figure(2); 

set (gcf, 'DefaultTextFontSize', fonts, 'DefaultAxesFontSize', fonts);
n_freq=100;f=freq(n_freq); k=k_real;
I_abs=abs(PHI0(n_freq,:));
semilogx(k, I_abs, 'linewidth',lines);
x_min=(k(2)-k(1));x_max=0.6*max(k);
x_min=.1;x_max=1;

y_max=10*floor(max(I_abs/5))+5;y_min=0;dy=0.06*y_max;
semilogx(real(k),I_abs,'b','linewidth',lines);grid;
axis( [x_min  x_max y_min y_max]); 
grid;
xlabel('Horizontal wave number – 1/m')
ylabel('Integrand');
text(x_min, y_max-dy, ['Freq.= ' num2str(f) ' Hz']);
text(x_min, y_max-3*dy, [' D =   ' num2str(D) ' m']);
text(x_min, y_max-4*dy, [' SD =   ' num2str(zs) ' m']);
text(x_min, y_max-5*dy, [' RD =   ' num2str(z) ' m']);

 
