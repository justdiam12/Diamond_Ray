
%%   Author: Jens M. Hovem
%%   Copyright 2011 Acoustic Research Center, NTNU
%%   Revised $  $Date: 2011/01/11

%%   Revised Date: 2012/01/21
%plot_time_signals
disp( 'loading the files from disk "t_real","signal" and "phone_positions"')
load signal;load phone_positions; load t_real;

fs=para.fs;
z_source=para.z_source; z_receiver=para.z_receiver;
phone_positions = newpara.range_phone;
newpara.range_receiver = newpara.range_phone;
N_phones = length(phone_positions);
phone_spacing=mean(phone_positions(2:end)-phone_positions(1:end-1));
rmax=max(phone_positions)+ 2*phone_spacing;
rmin=min(phone_positions)-phone_spacing;
max_signal=max(max(signal));
signal=signal/max_signal;
% if ~isnan(phone_spacing)
%     scale=phone_spacing/(1000*max_signal);
% else
%     scale=1;
% end;
% max_signal=max(max(signal));
% 
% if ~isnan(phone_spacing)
%     scale=phone_spacing/(1000*max_signal);
% else
%     scale=1;
% end;
fonts=14; lines=1;
% Time responses plotted in real time
figure(20);clf;
disp('Figure 20 plots the total time responses in real time');
scale=1;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

for jj=1:1: N_phones
    hold on;
    plot(t_real(jj,:),signal(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines);
end;
grid
xlabel('Time – s'); ylabel ('Range – km');
x_min=min(min(t_real));
x_max=max(max(t_real));
y_min=rmin/1000;
y_max=rmax/1000;
axis([x_min x_max y_min y_max]);
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])
box on

% Plot the "other arrivals" by setting other=1;
other =0;
figure(21); clf;
disp( 'Figure 21 plots time responses  in real time and sorted')

set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
scale=scale*7;
for jj=1:N_phones
    scale=1000;
    hold on;
    t_plot=t_real(jj,:);
    if other==1;plot(t_plot,signal_other(jj,:)*scale+bias(jj)/1000,'y','LineWidth',lines); end;%other
    plot(t_plot,signal_direct(jj,:)*scale+bias(jj)/1000,'k','LineWidth',lines); % direct signal
    plot(t_plot,signal_bottom(jj,:)*scale+bias(jj)/1000,'m','LineWidth',lines);% bottom
    plot(t_plot,signal_surface(jj,:)*scale+bias(jj)/1000,'g','LineWidth',lines);%surface
    plot(t_plot,signal_turn(jj,:)*scale+bias(jj)/1000,'r','LineWidth',lines);%refracted
    plot(t_plot,signal_sufbot(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines); %surface-bottom
    
end;
x_min=min(min(t_real));
x_max=max(max(t_real));
y_min=rmin/1000;
y_max=rmax/1000;
%axis([ x_min x_max y_min y_max])
text(0.75* x_max, 0.1*y_max, ['{c_{red}}=' num2str(fix(para.c_red)) ' m/s'])
if other==1;legend('Other', 'Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end
if other==0;legend('Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end ;
xlabel('Time – s'); ylabel ('Range – km');
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])
grid

figure(22); clf;
disp( 'Figure 22 plots time responses  in reduced time and sorted')

rr=1:para.Rmax;
c0=env.ssp(end,2);  c1=env.cp1; c_red=para.c_red;
theta=acos(c0/c1);

set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
for jj=1:N_phones
    hold on;
    t_plot=t_red(jj,:);
    if other==1; plot(t_plot,signal_other(jj,:)*scale+bias(jj)/1000,'y','LineWidth',lines); end;%other;
    plot(t_plot,signal_direct(jj,:)*scale+bias(jj)/1000,'k','LineWidth',lines); % direct signal
    plot(t_plot,signal_bottom(jj,:)*scale+bias(jj)/1000,'m','LineWidth',lines);% bottom
    plot(t_plot,signal_surface(jj,:)*scale+bias(jj)/1000,'g','LineWidth',lines);%surface
    plot(t_plot,signal_turn(jj,:)*scale+bias(jj)/1000,'r','LineWidth',lines);%refracted
    plot(t_plot,signal_sufbot(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines); %surface-bottom
    
end;
x_min=0; x_max=(N/fs)/8;
y_min=rmin/1000;
y_max=rmax/1000;
axis([ x_min x_max y_min y_max])
text(0.75* x_max, 0.1*y_max, ['{c_{red}}=' num2str(fix(para.c_red)) ' m/s'])
if other==1;legend('Other', 'Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end
if other==0;legend('Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end ;
xlabel('Reduced time – s'); ylabel ('Range – km');
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])
grid

figure(23); clf;
disp( 'Figure 23 plots time responses in reduced time and compensated for spherical spreading loss')

set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
for jj=1:N_phones
    hold on;
    factor=phone_positions(jj)/1000;
    sign_d=signal_direct*factor;
    sign_b=signal_bottom*factor;
    sign_s=signal_surface*factor;
    sign_t=signal_turn*factor;
    sign_sb=signal_sufbot*factor;
    sign_o=signal_other*factor;
    t_plot=t_red(jj,:);
    if other==1; plot(t_plot,sign_o(jj,:)*scale+bias(jj)/1000,'y','LineWidth',lines);end; %other
plot(t_plot,sign_d(jj,:)*scale+bias(jj)/1000,'k','LineWidth',lines); % direct signal;
plot(t_plot,sign_b(jj,:)*scale+bias(jj)/1000,'m','LineWidth',lines);% bottom
plot(t_plot,sign_s(jj,:)*scale+bias(jj)/1000,'g','LineWidth',lines);%surface
plot(t_plot,sign_t(jj,:)*scale+bias(jj)/1000,'r','LineWidth',lines);%refracted
plot(t_plot,sign_sb(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines); %surface-bottom

end;
axis([ x_min x_max y_min y_max])

box on

axis([ x_min x_max y_min y_max])
text(0.75* x_max, 0.1*y_max, ['{c_{red}}=' num2str(fix(para.c_red)) ' m/s'])
if other==1;legend('Other', 'Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end
if other==0;legend('Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end ;
xlabel('Reduced time – s'); ylabel ('Range – km');
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])


figure(24); clf;
disp( 'Figure 24 plots the absolute values of the time responses in reduced time and compensated for spherical spreading loss')

set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
for jj=1:N_phones
    hold on;
    factor=phone_positions(jj)/1000;
    
    sign_d=abs(signal_direct)*factor;
    sign_b=abs(signal_bottom)*factor;
    sign_s=abs(signal_surface)*factor;
    sign_t=abs(signal_turn)*factor;
    sign_sb=abs(signal_sufbot)*factor;
    sign_o=abs(signal_other)*factor;
    t_plot=t_red(jj,:);
    if other==1; plot(t_plot,sign_o(jj,:)*scale+bias(jj)/1000,'y','LineWidth',lines); end;%other
plot(t_plot,sign_d(jj,:)*scale+bias(jj)/1000,'k','LineWidth',lines); % direct signal;
plot(t_plot,sign_b(jj,:)*scale+bias(jj)/1000,'m','LineWidth',lines);% bottom
plot(t_plot,sign_s(jj,:)*scale+bias(jj)/1000,'g','LineWidth',lines);%surface
plot(t_plot,sign_t(jj,:)*scale+bias(jj)/1000,'r','LineWidth',lines);%refracted
plot(t_plot,sign_sb(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines); %surface-bottom


end;

axis([ x_min x_max y_min y_max])

text(0.75* x_max, 0.1*y_max, ['{c_{red}}=' num2str(fix(para.c_red)) ' m/s'])
if other==1;legend('Other', 'Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end
if other==0;legend('Direct','Bottom', 'Surface','Refracted','Surface-Bottom');end ;
xlabel('Reduced time – s'); ylabel ('Range – km')
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])



figure(25); clf;
disp( 'Figure 25 plots the total time responses in reduced time ')

%scale=scale/3;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
for jj=1:N_phones
    hold on;
    factor=phone_positions(jj)/1000;
    sign_d=abs(signal_direct)*factor;
    sign_b=abs(signal_bottom)*factor;
    sign_s=abs(signal_surface)*factor;
    sign_t=abs(signal_turn)*factor;
    sign_sb=abs(signal_sufbot)*factor;
    t_plot=t_red(jj,:);
    plot(1000*t_plot,signal(jj,:)*scale+bias(jj)/1000,'b','LineWidth',lines);
    
end;

x_min=0; x_max=(N/fs*.5)*1000;
y_min=rmin/1000;
y_max=rmax/1000;
%axis([ x_min x_max y_min y_max])

text(0.75* x_max, 0.1*y_max, ['{c_{red}}=' num2str(fix(para.c_red)) ' m/s'])
legend('Total response');
xlabel('Reduced time – s'); ylabel ('Range – km')
title([env.title ':   Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])

