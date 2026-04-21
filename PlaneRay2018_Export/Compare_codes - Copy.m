
close all;
%clear
lines=1.0;fonts=16;
load TL_plr; load r_plr;
load para; load env;
freq=para.frequency;

%freq=[ 25 50 75 100 125 150 175 200];
M=length(freq);
zs=para.z_source; zr=para.z_receiver;
Type= input('Give the input number')
if Type==1;
    title_text='Fluid bottom';
    load RangeOast1;load TLOast1;
    roast=RangeOast1;
    TLoast=TLOast1;
end

if Type==3;
    title_text='Layered bottom';
    load RangeOast3;load TLOast3;
    roast=RangeOast3;
    TLoast=TLOast3;
end


for m=1:M;
    figure(m);
    fr=para.frequency(m);
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    r=roast(:,1,m);tl=TLoast(:, 1,m);
    tl=tl(r>0); r=r(r>0);
    plot(r, tl, 'b','linewidth',lines);
    
    axis([0 5 20 120]);grid
    axis ij;legend( 'OASES' )
    text(0, 22,['  Freq = ' num2str(fr) ' Hz' ])
    text(0, 27,['  SD =  ' num2str(zs) ' m' ]);
    text(0, 32,['  SR =  ' num2str(zr) ' m' ]);
    title([  title_text])
end
legend( 'OASES' )
ylabel('Transmission loss - dB')
xlabel('Range - km')
title([ title_text])
text(0, 22,['  Freq = ' num2str(fr) ' Hz' ])
text(0, 27,['  SD =  ' num2str(zs) ' m' ])
text(0, 32,['  SR =  ' num2str(zr) ' m' ]);





