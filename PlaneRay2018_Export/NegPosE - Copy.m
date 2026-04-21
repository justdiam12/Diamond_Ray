%NegPos
%title_tekst='Negative&Positive gradient'
clear
        del_z=.1; wd=200;
        z_source=25; z_receiver=wd;
        z_break=50;
        z=0:del_z:wd;
        g1=-0.02; g2=0.1; c0=1480;
                z1=1:z_break;
        z2=z_break:wd;
        
        c1=c0+g1*z1;
        c2=c0+g1*z_break+g2*(z2- z_break);
        c=[c1 c2]; z=[z1 z2];
        
       % c( z>z_break)=c1+g2*z(z> z_break)+ g1*z_break;
        
        figure(88) ;clf;
        plot(c,z,'-');axis ij;
        axis([ 1470 1500 0 200])