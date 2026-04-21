
%   Author: Jens M. Hovem 
%   Revised $  $Date: 2014/06/15 
%   Other start pulses may be added


function [source_signal, fs, t_start] = getsourcesignal(source_type,para)
N=para.nfft;n=1:N;
fs=para.fs;
 source_signal=zeros(1,N);

if source_type==1; %Specifications for the Ricker source pulse
    f_max=fs/10;
    T_D=sqrt(20)/(pi*f_max);% Time separation between the two negative peaks
    t_start=1*T_D;
    source=Ricker(f_max,t_start,fs,N);
    source_signal=source;
end

if source_type==2; %Specifications of a Gaussian  pulse
    T=.01; %Pulse duration i sec;
    sigma_t=T/2;
    t=(0:N-1)/fs;
    arg=((t-T).^2/sigma_t^2);
        amp=1;
    pulse=amp*exp(-arg);
    source_signal=pulse;
end   
if source_type==3; %Specifications for a "delta" pulse
    t_start=0.01;
    N_start=fix(t_start*fs);
    source_signal(N_start)=1;
end







            
            