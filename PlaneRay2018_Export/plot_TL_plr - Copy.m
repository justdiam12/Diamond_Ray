%plot_Tl_plr

close all;clear;
lines=1.5;fonts=16;
select_freq=[ 25 50 75  100 125 150 175  200];
q=[ 1 2 3 4 5 6 7 8];
M=length(select_freq);
load c:\Planeray2018\TL_plr;
load c:\Planeray2018\r_plr;
load c:\Planeray2018\env;
load c:\Planeray2018\para;




[R, Q]=size(TL_plr);
for q=1:1:Q
    SD=para.z_source; RD=para.z_receiver;freq=para.frequency;
    WD= max(env.bathy(2,1));
    figure(9+q);
        fonts=12;lines=1.5;
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    maxy=120;miny=20;maxx=para.Rmax/1000; minx=0;  
    plot(r_plr/1000,TL_plr(:,q),'g', 'linewidth',lines);
    axis ([ 0 max(r_plr)/1000 20 120]);axis ij
    xlabel ('Range - km'); ylabel( 'Transmission loss - dB ');
    title_text=env.title;
   text(0, miny+5,[ title_text   ': Freq = ' num2str(freq(q)) ' Hz']);
   text(maxx*0.7, miny+5,[' WD = ' num2str(WD) ' m']);
   text(maxx*0.7, miny+10,[' SD = ' num2str(SD) ' m'])
   text(maxx*0.7, miny+15,[' RD = ' num2str(RD) ' m'])
   grid
   
   
end