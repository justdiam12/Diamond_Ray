%planeray_depth

%%   Author: Jens M. Hovem
%%   Copyright 2011 Acoustic Research Center, NTNU
%%   Revised:  2011/01/11 .
%5   This program contains a receiver-depth loop and displays received time signals
%%   as function of receiver depth for a specific range.
%%   The program also calulates the peak pressure level(PPL) and sound exposure
%%   level(SEL) for each of the timetraces and plos the rersult as function of depth
%%   The source level is given by the parameter SL (dB rel. 1 microPa, at 1 m.

close all;clear all;

Example=input( 'Example ? ');
rawinput = planerayinput(Example);
[env, para] = initpara(rawinput);
save env; save para;
wd=max(env.bathy(2,:));
SL=200;

range=input('Specify the ranges  you want to study......');

[temp1,range_index] = min(abs(range- para.range_receiver));
para.range_phone=range;
receiver_depth=input('Specify the depths you want to study......');
receiver_depth=receiver_depth(receiver_depth<wd);
N=para.nfft; fs=para.fs;
n=0:N-1; time=n/fs;
para.timeplot=0; para.timeplot=0;para.envplot=0; %Turn off the plotting of time plots to save time

[source_signal, t_start] = getsourcesignal(1,para);
Q=length(receiver_depth);

ind=zeros(1,Q);legendstr=cell(1,Q);
para.envplot=1;
para.quiet=1;
for q=1:Q
    disp(['Receiver depth = ' num2str(receiver_depth(q)) ' m']);
    para.z_receiver=receiver_depth(q);
    r=para.range_receiver;
    [HIST, COUNT, SUFBOT, Rays] = tracerays(env, para);
    [eigenangle, HIST, COUNT, SUFBOT, Rays, para]=...........
        sortrays(HIST, COUNT, SUFBOT, Rays, env, para);
    [transfer_function,trans,c_red,eigen]  = transfunc(HIST, COUNT, SUFBOT, env, para,1);
    trans_loss = -20*log10(eps^2+abs(trans.total));
    TL_depth(q,:)= trans_loss(:,1);
    depth=receiver_depth(q);
    legendstr{q}=['RD = ',num2str(depth),' m'];
    RD(q,:)=depth;
    EIG=HIST;
    para.z_receiver=receiver_depth(q);
     [pulse, h, eigen]=timeresponse(EIG, COUNT, SUFBOT, env, para, source_signal, t_start);
    signal(q,:)=pulse(1:para.nfft);
end;

figure(70); clf;
fonts=16; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
plot(r/1000, TL_depth,'linewidth',lines);
grid; legend(legendstr);
xlabel('Range - km');ylabel(' Transmission loss -dB');
axis ij; axis([0 5 20 120 ])
  


figure(71); clf; 
fonts=16 ;lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

plot(SL-TL_depth(:,range_index), receiver_depth,'linewidth', lines );grid
xlabel('Sound level - dB')
ylabel(' Receiver depth -m ')
xmin=100; xmax=160; ymin=0;ymax=max(receiver_depth);
axis ij; axis([ xmin xmax ymin ymax])
text( xmin, ymin +15,['  SL= ' num2str(SL) , ' dB'])
text( xmin, ymin +30,['  SD= ' num2str( para.z_source) , ' m'])
text( xmin, ymin +45,['  Range= ' num2str( range/1000) , ' km'])
text( xmin, ymin +60,['  Freq ' num2str(para.frequency) , ' Hz'])
%text( xmin, ymin +40,['  SL= ' num2str(SL) , 'dB'])


figure(72); clf;
fonts=16 ;lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

plot(TL_depth(:,range_index), receiver_depth,'linewidth', lines );grid
xlabel('Transmission loss - dB')
ylabel(' Receiver depth -m ')
axis ij; axis([ 20 120 0 max(receiver_depth)])
title([env.title ':   Sd = ' num2str(para.z_source) ' m, ' ' Rd= = ' num2str(RD(q)) ' m' ]);

figure(73); clf;
fonts=16 ;lines=2;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts)

 scale=10000;
x_min=0; x_max=max(time)*1000;

for q=1:Q; plot(time,scale*signal(q,:)+10*q,'linewidth',lines);hold on; end
axis ij;axis([ 0 0.25 0 250])
xlabel( 'Reduced time - ms'); ylabel( 'Depth - m')
title(['Time responses: ' env.title   ' Sd = ' num2str(para.z_source) ' m ' ' Range = ' num2str(range/1000) ' km' ]);

figure(74);
fonts=16 ;lines=2;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts)
q=3
plot(r/1000,TL_depth(q,:), 'linewidth', lines);grid
axis ij;axis([ 0 5 20 120])
xlabel('Range'); ylabel('Transmission loss - dB')
title([env.title ':   Sd = ' num2str(para.z_source) ' m, ' ' Rd= = ' num2str(RD(q)) ' m' ]);


figure(75); clf
fonts=16 ;lines=2;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts)
TL= TL_depth;
ma=max(max(TL)); mi=min(min(TL));
ma =5; mi=100;
mst=(ma-mi)/80;
cvec=[mi:mst:ma];
plotssp=0;
colormap(flipud(jet))
contourf(r/1000,RD,TL,cvec);axis ij;colorbar; shading flat;

hold on;
plot_bathy(env,para, plotssp)
xlabel('Range - km '); ylabel('Depth - m')
title([env.title ':   Sd = ' num2str(para.z_source) ' m, ' ' Freq  = ' num2str(para.frequency) ' Hz' ]);


