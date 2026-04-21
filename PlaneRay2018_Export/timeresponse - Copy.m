%%   Authors: Jens M. Hovem and Shefeng Yan
%% July 2013
 
function [signal_conv, h, eigen] = timeresponse(HIST, COUNT,SUFBOT, env, para, source_signal, ~)
phone_positions=para.range_phone;
newpara=para;
if length(phone_positions)>1000;
    disp( 'Too many phone positions, reduced to 1000')
    range_phone= linspace(min(phone_positions),max(phone_positions),1000);
    newpara.range_phone=range_phone;
end;
fs=para.fs;
%z_source=para.z_source; z_receiver=para.z_receiver;
phone_positions = newpara.range_phone;
newpara.range_receiver = newpara.range_phone;
N_phones = length(phone_positions);
%phone_spacing=mean(phone_positions(2:end)-phone_positions(1:end-1));
%f_max=fs/10;
%t_start=5/f_max;

N= para.nfft;
W=fft(source_signal,N);
source_spectrum=W(1:N/2);
t=(0:N-1)/fs;
f=(0:N/2-1)*fs/N;
  newpara.frequency=f;
    %The Source signal is stored for future use
    save source_signal source_signal;



[transfer_function,transfer,c_red,eigen]  = transfunc(HIST, COUNT, SUFBOT, env, newpara,2);
save newpara
save t

signal=zeros(N_phones,para.nfft);
signal_direct=signal;signal_turn=signal;signal_surface=signal;
signal_bottom=signal;signal_sufbot=signal;signal_other=signal;
bias=signal;factor=zeros(N_phones); t_red=signal;t_real=signal;h=signal;
signal_filter=signal;
%signal_conv=signal;
for jj=1:N_phones
    receiver_range=phone_positions(jj);
    
    
    if (max(eigen.timedelay(jj,:),[],2)-receiver_range/c_red) >= N/fs
        warning(['N/fs too small when receiver range is ',num2str(receiver_range),' m']);
    end;
    
    % Convert to time domain
    signal_spectrum=source_spectrum.*transfer_function(jj,:);
    signal_spectrum_direct=source_spectrum.*transfer.direct(jj,:);
    signal_spectrum_turn=source_spectrum.*transfer.turn(jj,:);
    signal_spectrum_surface=source_spectrum.*transfer.surface(jj,:);
    signal_spectrum_bottom=source_spectrum.*transfer.bottom(jj,:);
    signal_spectrum_sufbot=source_spectrum.*transfer.sufbot(jj,:);
    signal_spectrum_other=source_spectrum.*transfer.other(jj,:);
    signal(jj,:)=2*real(ifft(signal_spectrum,N)).';
    signal_direct(jj,:)=2*real(ifft(signal_spectrum_direct,N)).';
    signal_turn(jj,:)=2*real(ifft(signal_spectrum_turn,N)).';
    signal_surface(jj,:)=2*real(ifft(signal_spectrum_surface,N)).';
    signal_bottom(jj,:)=2*real(ifft(signal_spectrum_bottom,N)).';
    signal_sufbot(jj,:)=2*real(ifft(signal_spectrum_sufbot,N)).';
    signal_other(jj,:)=2*real(ifft(signal_spectrum_other,N)).';
    bias(jj)=phone_positions(jj);
    factor(jj)=phone_positions(jj)/1000;
    
    
    t_red(jj,:)=t;
    t_real(jj,:)=t+(receiver_range-para.range_source)/c_red;
    H=[transfer_function(jj,:),0,fliplr(conj(transfer_function(jj,2:end)))];
    h(jj,:)=real(ifft(H));
    signal_filter(jj,:)=filter(h(jj,:),1,source_signal);
   signal_conv(jj,:)=conv(h(jj,:),source_signal);
end; %jj over receiver positions
save signal signal
save phone_positions phone_positions
save t_real;
save t_red;

if~para.quiet
 disp( 'The results are stored on disk as "signal" and "phone_positions"')
end;





