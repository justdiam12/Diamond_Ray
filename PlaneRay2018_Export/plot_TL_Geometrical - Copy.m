% plot_TL_Geometrical

fonts=16; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
wd=200;
r=1:100000;
r0=100;
d1=100; D=d1-0.01*r;
[TL] = TL_geometrical_depth(r, r0, D);
semilogx(r/1000, TL, 'linewidth',lines)
xlabel( 'Range km'); ylabel( 'Transmision loss dB')