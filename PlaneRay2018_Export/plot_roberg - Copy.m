%plot_roberg
close all
figure(1)
fonts=20 ;lines=2;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

%load ssp_roberg;load z_roberg;

%Q=length(para.frequency);
Q=12;
ind=zeros(1,Q);legendstr=cell(1,Q);
for q=1:Q; 
        legendstr{q}=['Month = ',num2str(q)];    
end;

plot( ssp_roberg', z_roberg', 'linewidth', lines);axis ij;grid;
xlabel('Sound speed – m/s ')
ylabel('Depth – m/s')
legend(legendstr)
axis( [1460 1520 0 500])
title('Seasonal variations - Trondheim fjord')

figure(2)

fonts=20 ;lines=2;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

%load ssp_roberg;load z_roberg;
plot( ssp_roberg(2,:)', z_roberg(2,:), ssp_roberg(4,:)', z_roberg(4,:),ssp_roberg(6,:)', z_roberg(6,:),.............
    ssp_roberg(8,:)', z_roberg(8,:)',  ssp_roberg(10,:), z_roberg(10,:),'linewidth', lines)
axis ij;grid;
xlabel('Sound speed – m/s ')
ylabel('Depth – m/s')

%axis( [1470 1500 0 200])
legend(' February', ' April', 'June',' August', ' October', 'Location', 'SouthEast')
title('Seasonal variations - Trondheim fjord')

 