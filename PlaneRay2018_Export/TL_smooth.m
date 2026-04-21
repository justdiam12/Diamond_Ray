%TL_smooth
clear



figure
lines =2;fonts=18;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
load para; load env;load TL_plr; load r_plr
SD=para.z_source; RD=para.z_receiver;
wd=env.bathy(2,1);
titlet_tekst=env.title;

TL_geo= TL_geometrical(r_plr, wd);
N=10;n=1:N;
A=1;B=(1/N)*ones(1,N);
S=10.^(-TL_plr/10);
Y = filter(B,A,S);

TL_av=-10*log10(Y);
semilogx(r_plr,TL_av,r_plr, TL_geo,'--', 'linewidth',lines);axis ij;
xlabel( 'Range - m'); ylabel('Transmission loss - dB')

%title(title_tekst)
text (101,5,[ 'WD = ' num2str(wd) ' m'])
text (101,15,[ 'SD = ' num2str(SD) ' m '])
text (101,25,[ 'RD = ' num2str(RD) ' m'])
axis([100 10000 0 100 ])   
title(env.title)
Q=length(para.frequency);
ind=zeros(1,Q);legendstr=cell(1,Q);
for q=1:Q; 
        legendstr{q}=['Freq = ',num2str(para.frequency(q)),' Hz'];    
end;

legend( [ legendstr 'Geometrical']) ;grid