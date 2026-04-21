function [caustic_range , caustic_amp ] = FindCaustics(HIST, eigen,para, env, plotting )

%FindCaustics
%%  Date: 2015/03/25
%% This program is derived from the program planerayfield used planeray
%% Input; HIST,COUNT, para,env; most recent PlaneRay results

%% Output  caustics=[caustic_range' caustic_angle' caustic_amp'];
%% figure(plots): Ray history and initial angles of the peaks
%% figure(plots+1)Ray amplitudes of peaks as function of range

if nargin<5;plotting=0;end;
RANGE=HIST.RANGE;
THETA=HIST.THETA*180/pi;
amplitude=eigen. rayamp;
theta=eigen.startangle;
theta_max=max(max(abs(theta)))*0.5;
range=eigen.range;
[~, K]= size(amplitude);
if K>1;
    Z=zeros(K);
    caustic_range=Z;    caustic_angle=Z;    caustic_amp=Z;
    q=0;
    for k=1:K;
        amp= abs(amplitude(:,k));
        [ray_amp(k), m]=max(amp);
        temp_angle=theta(m,k);
        temp_range=range(m);

        if temp_range>0&& abs(temp_angle)<theta_max; q=q+1;
            caustic_range(q)=temp_range;
            caustic_angle(q)=temp_angle;
            caustic_amp(q)=ray_amp(k);
        end
        caustics=[caustic_range' caustic_amp'];
    end;
    
    for k=1:q;
        %disp([' N_caustic = ' num2str( k) ',  Range '  num2str(fix(caustic_range(k))) ' m, ' ' Initial angle  ' num2str(caustic_angle(k)) ',  deg']);
         disp([' N_caustic = ' num2str( k) ',  Range '  num2str(caustic_range(k)) ' m, ' ' Initial angle  ' num2str(caustic_angle(k)) ',  deg']);
    end
    if plotting>0;
        figure(plotting); clf; fonts=20; lines=2;
        set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
        
        plot( caustic_range/1000,caustic_angle,'ro',RANGE/1000, THETA,'linewidth',lines);axis ij
        title([ env.title  ':  Sd=' num2str(para.z_source) ' m, ','Rd=' num2str(para.z_receiver) ' m ']);
        start_angle=para.start_angles;
        x_min=eps;x_max= para.Rmax/1000; y_min=min(start_angle-5); y_max=max(start_angle+5);
        axis([x_min x_max y_min y_max]); grid
        xlabel('Range – km'); ylabel('Initial angle  – deg')
       else
        disp(' No caustics found');
    end
end;
figure(plotting+1);  
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
amplitude=amplitude/ max(max(amplitude));

plot(range/1000,amplitude,'-o','linewidth', lines);
title([ env.title  ':  Sd=' num2str(para.z_source) ' m, ','Rd=' num2str(para.z_receiver) ' m ']);

x_min=0;x_max= para.Rmax/1000; y_min=0; y_max=1.2;
axis([x_min x_max y_min y_max]); grid
xlabel('Range – km'); ylabel('Ray amplitude')



end





