function ploteigstructure(Rays, EIG, COUNT, env, para)

%   Produces plots of range to receiver array as function of initial angle,
%   and geometrical spreding loss as function of range for each ray class
%   Authors: Shefeng Yan and Jens M. Hovem
%   
%   $Revised   $ Date: 2013/28/07


eigen = planerayfield(EIG, COUNT, env,para);
rayamp=eigen.rayamp;
count=eigen.count;

Rmax=para.Rmax;
eig_range = EIG.RANGE;
eig_theta = EIG.THETA;

eig_angle=EIG.ANGLE;
eig_time=EIG.TIME;
eig_theta=eig_theta*180/pi;%Convert to degrees
%eig_angle=eig_angle*180/pi;%Convert to degrees
attplot=0;
if any(any(eig_theta==0)==1)
    eig_theta(eig_theta==0)=nan;
    attplot=1;
end;

transloss = -20*log10(rayamp+eps);
transloss(isinf(transloss))=200;
save transloss transloss;
save eig_range eig_range;





figure(30); clf
if ~para.quiet; disp(' Figure 30 plots range vs initial source angle'); end;
fonts=24; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);


plot(eig_range/1000,eig_theta, '.',eig_range/1000,eig_theta, '- ','linewidth',lines );
grid
title([ env.title  ':  Sd=' num2str(para.z_source) ' m, ','Rd=' num2str(para.z_receiver) ' m ']);
start_angle=para.start_angles;
x_min=0;x_max= Rmax/1000; y_min=min(start_angle-5); y_max=max(start_angle+5);
axis([x_min x_max y_min y_max]);
xlabel('Range – km');
ylabel('Initial angle  – deg')
set(gca,'ydir','reverse');


figure(31); clf
if ~para.quiet;disp('Figure 31 plots travel time as function of initial source angle '); end
fonts=24; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

plot(eig_time,eig_theta, '.',eig_time,eig_theta,'-',  'linewidth',lines);

grid
title([ env.title  ':  Sd=' num2str(para.z_source) ' m, ','Rd=' num2str(para.z_receiver) ' m ']);

x_min=0;x_max= max(max(eig_time+0.2)); y_min=min(start_angle-5); y_max=max(start_angle+5);
axis([x_min x_max y_min y_max]);

xlabel('Travel time – s');
ylabel('Initial angle  – deg')
set(gca,'ydir','reverse');




figure (32); clf;
if ~para.quiet;disp('Figure 32 plots geometrical transmission loss for the various ray classes');end
fonts=16; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);%plots Range vs initilal angle at reciver
plot(para.range_receiver/1000,transloss,'linewidth',lines); grid;axis ij;
title([ env.title  ':  Sd=' num2str(para.z_source) ' m, ','Rd=' num2str(para.z_receiver) ' m ']);
axis([ 0 para.Rmax/1000 0 120]);
xlabel('Range – km');
ylabel('Geom. tran. loss – dB')




