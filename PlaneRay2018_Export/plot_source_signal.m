function [source_spectrum, freq, time ]= plot_source_signal(source_signal,para, figure_no)
%   Plot the shape and the spectrum of the source siganl;
%   Author: Jens M. Hovem
%   $Revision: 3.0 $  $Date: 2014/06/15$
load para
fs=para.fs;
N=length(source_signal);

W=fft(source_signal,N);
source_spectrum=W(1:N/2);
time=(0:N-1)/fs;
freq=(0:N/2-1)*fs/N;
if figure_no>0;
    figure(figure_no);clf;fonts=18; lines=2;
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    subplot(2,1,1);
    plot(time*1000,source_signal,'k','linewidth',lines)
    xlabel( 'Time 	– ms','Fontsize',fonts)
    ylabel('Amplitude','Fontsize',fonts);
    axis([ 0 max(time)*1000/5 -1.5  1.5])
    subplot(2,1,2)
    plot(freq(1:N/2),abs(source_spectrum),'k','linewidth',lines);
    xlabel(' Frquency – Hz','Fontsize',fonts)
    ylabel(' Spectrum amplitude','Fontsize',fonts)
end
