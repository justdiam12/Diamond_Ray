function SEL =compute_SEL(trans,source_signal,env,para);
% Calculate the Sound Exposure Level by integration in the frequency domain
% The calculations are based on the store results of the last run of PlaneRay
%   Author: Jens M. Hovem 
%   Copyright 2011 Acoustic Research Center, NTNU
%   $Revision: 1.0 $  $Date: 2011/05/10 
   

N=para.nfft;
phone_positions=para.range_phone;
N_phones=length(phone_positions); 
source_spectrum=fft(source_signal, N);
pulse_spectrum=zeros( N_phones,N/2);
%SEL_freq=zeros(N_phones);
for n=1:N_phones
    
    transfer_function=trans.total(n,:);
    pulse_spectrum(n,:)= transfer_function.*source_spectrum(1:N/2);% One_sided spectrum
    freq_spec=pulse_spectrum(n,:);
    freq_spec(1)=0.5*freq_spec(1);
    freq_spec(N/2)=0.5*freq_spec(N/2);
    power_freq =4*sum(real(freq_spec).^2/N);
    SEL_freq(n)=10*log10((power_freq)/N);
end
SL=220;
SEL=SEL_freq;
semilogx(phone_positions, SL+SEL)