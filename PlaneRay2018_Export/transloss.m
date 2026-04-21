function [transloss_dB,range] =transloss(HIST, COUNT, SUFBOT, env, para, plots)
%  Calculates and plot the transmission loss for the frequencies specified 
%  in para.frequency 
%   Authors: Jens M. Hovem and Shefeng Yan
%   July 2013
%plots>1; linear range scale; plots<0 logaritmic range scale; 
nfft=para.nfft;
fs=para.fs;
carrier_frequency=para.carrier_frequency;
f_select=para.frequency;
z_source=para.z_source; z_receiver=para.z_receiver;

range=para.range_receiver;

conplot=para.conplot;
newpara=para;
if conplot && nfft~=0 && fs~=0
    newpara.frequency = (0:nfft/2-1)*fs/nfft;
end;
[transfer_function,trans,c_red,eigen]  = transfunc(HIST, COUNT, SUFBOT, env, para,1);


trans_loss = -20*log10(eps^2+abs(trans.total));
Q=length(f_select);
ind=zeros(1,Q);legendstr=cell(1,Q);
for q=1:Q
    freq=f_select(q);
    [temp1,ind(q)] = min(abs(freq-newpara.frequency));
    if carrier_frequency==0;
        legendstr{q}=['Freq = ',num2str(newpara.frequency(ind(q))),' Hz'];
    else
        legendstr{q}=['Freq = ',num2str((newpara.frequency(ind(q))+carrier_frequency)/1000),' kHz'];
    end;
    
end;

TLdB=trans_loss(:,ind);
disp('TL_plr and r_plr are stored on the  disk for later use')
TL_plr=trans_loss;
r_plr=range;
transloss_dB=TLdB;
save TL_plr TL_plr;
save r_plr r_plr;



if para.tlplot==1  %Start plotting section
     figure(10); clf %
    %disp('Figure 10 plots tranmission loss on a linear range scale')
    fonts=16; lines=1.5;
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    x_min=min(range/1000); x_max=para.Rmax/1000;
    dB_min=20; dB_max=120;
    plot(range/1000,TLdB,'linewidth',lines);  grid on
    set(gca,'ydir','reverse');
    axis([x_min x_max dB_min dB_max])
    legend(legendstr);
    xlabel('Range – km');
    ylabel('Transmission loss – dB');
    title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m ']);
    
    figure(11); clf %
    %disp('Figure 10 plots tranmission loss on a logaritmic  range scale')
  
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    x_min=min(range/1000); x_max=para.Rmax/1000;
    dB_min=20; dB_max=120;
    
    semilogx(range/1000,TLdB,'linewidth',lines); grid on
    set(gca,'ydir','reverse');
    axis([x_min x_max dB_min dB_max])
    legend(legendstr);
    xlabel('Range – km');
    ylabel('Transmission loss – dB');
    title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m ']);
    
    


 if conplot;
    disp('Figure 11 Contoured transission loss as function og frequency and range');
    figure(11);
    fonts=12;
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    R=range(1:length(range));
    V=[20 30 40 50 60 70 80 90 100];
    
    contourf(newpara.frequency,R/1000,trans_loss,V);
    
    axis([0 fs/2 min(R)/1000-0.100, max(R)/1000+0.100])
    xlabel('Frequency – Hz')
    ylabel('Range – km');
     title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m ']);
    colormap(flipud(jet))
    colorbar
 end;
end % End plotting section

