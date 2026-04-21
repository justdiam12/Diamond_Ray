function p= ricker(f_0,t_0,f_s,N);
% Generates a Ricker pulse with
% f_s Sampling frequency in Hz
% f_0 frequency of maximum amplitude in Hz
% t_0= delay for time of the main peak of the Ricker wavelet. 
% N block length , normally a power of two

f_m= f_0;
del_t=1/f_s;
n=1:N;
t=n*del_t;
f=(n/N)*f_s;
T_c=t_0;
w1=1-2*pi^2*f_m.^2*(t-T_c).^2;
w2=exp(-(pi*f_m*(t-T_c)).^2);
p=w1.*w2;

