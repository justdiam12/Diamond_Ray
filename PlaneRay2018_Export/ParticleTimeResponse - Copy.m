
%% The program ParticleTimeResponses is created by Jens M. Hovem
%% Date: 2016/09/08
%% Program computes the timeresponses of horizontal and vertical components of
%% particle velocity and the sound pressure as function of range
% Initial settinge are results of last run of planeray
%% The program is a derived from the function tranfunc.m and can, in the future, be integrated as an option

close all
new =1;
% Initial settinge based on results of last run of planeray
if new==1
    load SUFBOT; load HIST; load COUNT; load env; load para;
end;

BOTTOM_ANGLE=SUFBOT.BOTTOM_ANGLE;
BOTTOM_RANGE=SUFBOT.BOTTOM_RANGE;
SURFACE_ANGLE=SUFBOT.SURFACE_ANGLE;
SURFACE_RANGE=SUFBOT.SURFACE_RANGE;
origin=1;
eigen = planerayfield(HIST, COUNT, env,para);
threshold=0;
%[ eigresult] = validate( newpara, env, SUFBOT,eigen, threshold);

startangle=eigen.startangle;targetangle=eigen.targetangle;
timedelay=eigen.timedelay; rayamp=eigen.rayamp;amplitude=eigresult.amplitude;
distance=eigen.distance;count=eigen.count;
phone_positions=para.range_phone;
range=para.range_phone/1000;
N_phones=length(phone_positions);

R_bathy = env.bathy(1,:);Z_bathy = env.bathy(2,:);
z=env.ssp(:,1);c=env.ssp(:,2);del_z=z(2)-z(1);
start_angle=para.start_angles;
bottom_depth=Z_bathy(1);c_depth=c(fix(bottom_depth/del_z)-1);
c_red=para.c_red;konst=8.686*2*pi;
nfft=para.nfft;


sigma_surface = env.sigma_surface; %Sea surface roughness
sigma_bottom = env.sigma_bottom; %Bottom roughness
dB_lambda_p1=env.ap1;  alpha_p1=dB_lambda_p1/konst;
dB_lambda_p2=env.ap2;  alpha_p2=dB_lambda_p2/konst;
dB_lambda_s2=env.as2;  alpha_s2=dB_lambda_s2/konst;
ccp1=env.cp1*(1+1i*alpha_p1);ccp2=env.cp2*(1+1i*alpha_p2);
ccs2=env.cs2*(1+1i*alpha_s2);rrho1=env.rho1;
rrho2=env.rho2;RR=env.RR;sedthick=env.laythick;

%% End of initial settings
% Plotting the source pulse and spectrum if source_plot >0.
source_plot=0;
if source_plot>0;
    [source_signal, fs, t_start] = getsourcesignal(1,para);
    [source_spectrum, freq, time ]=plot_source_signal(source_signal,para, source_plot);
end;
N=para.nfft; fs=para.fs;

freq=(0:N/2-1)*fs/N; f=freq;
nfreq=length(f);
omega=2*pi*f;
t=(0:N-1)/fs;
transfer_function = zeros(N_phones,nfreq);
transfer_function_sin = zeros(N_phones,nfreq);
transfer_function_cos = zeros(N_phones,nfreq);
cp1=env.cp1;cp2=env.cp2;cs2=env.cp2;
ap1=env.ap1;ap2=env.ap2;as2=env.as2;
rho1=env.rho1;rho2=env.rho2;
T=env.temperature;S=env.salinity;d=env.laythick;
ph=8;% fixed ph value)
D=mean([ para.z_receiver para.z_source]);
c0=1500; rho0=1000;
rhoc=c0*rho0;

[temp, a, a1, a2, a3] = frangar(T,S,D,ph,(f)/1000);
alpha=a/1000; lambda=1500./para.frequency;
for jj=1:N_phones
    [n rnr]=min(abs(para.range_receiver-para.range_phone(jj)));
    % [ targetangle(rnr,:)' targetangle(rnr,:)' rayamp(rnr,:)']
    
    receiver_range=phone_positions(jj);
    for cnt=1:size(startangle,2)
        current_angle=startangle(rnr,cnt);
        current_targetangle=targetangle(rnr,cnt);
        theta=current_angle*pi/180;
        current_angle;
        if ~isnan(current_angle)
            categ=count(:,cnt);
            N_bottom  = categ(1); N_surface = categ(2);
            N_above   = categ(3); N_below   = categ(4);
            if N_bottom>para.bothit,break;end;
            surface_ref=1; bottom_ref=1; aboveturn_ref=1; belowturn_ref=1;
            amp=rayamp(rnr,cnt);
            attenuation=10.^(-distance(rnr,cnt)*alpha/20);
            delay=timedelay(rnr,cnt)-(receiver_range- para. range_source)/c_red;
            E=exp(-1i*omega*delay);
            if N_surface>0
                theta_surface = interp1q(start_angle,SURFACE_ANGLE,current_angle);
                if sigma_surface~=0
                    for pp=1:N_surface
                        theta=theta_surface(pp);
                        R_s = R_coh( theta,sigma_surface, lambda );
                        surface_ref = surface_ref.*R_s;
                    end
                end;
                surface_ref = surface_ref*(-1)^N_surface;
                
            end;
            
            if N_bottom>0
                theta_bottom = interp1q(start_angle,BOTTOM_ANGLE,current_angle);
                range_bottom = interp1q(start_angle,BOTTOM_RANGE,current_angle);
                
                if length(R_bathy)<=2 && length(RR)<=2 && ~any(diff(theta_bottom(1:N_bottom)))
                    theta=theta_bottom(1);
                    if abs(theta)<pi;  R_b = R_bottom(theta,f,c_depth,cp1,cp2,cs2,rho1,rho2,d);end;
                    R_rough= R_coh( theta,sigma_bottom, lambda );
                    bottom_ref = (R_rough.*R_b).^(N_bottom);
                   % bottom_ref = (R_b).^(N_bottom);
                else
                    for pp=1:N_bottom
                        bottom_range=range_bottom(pp);
                        theta=theta_bottom(pp);
                        if length(R_bathy)>2 % not flat bottom
                            bottom_depth=Z_bathy(find(bottom_range-R_bathy>=0,1,'last'));
                            nn=(bottom_depth/del_z);
                            if  ~isempty(nn);   c_depth=c(fix(nn));end
                        end;
                        
                        if abs(theta)<pi; R_b = R_bottom(theta,f,c_depth,cp1,cp2,cs2,rho1,rho2,d);end
                        R_rough= R_coh( theta,sigma_bottom, lambda );
                        bottom_ref = bottom_ref.*R_b.*R_rough;
                    end
                end;
            end;
            if N_above>0
                aboveturn_ref=(-1i)^N_above;
            end;
            if N_below>0
                belowturn_ref=(-1i)^N_below;
            end;
            trfunc = E.*amp.*surface_ref.*bottom_ref.*aboveturn_ref.*belowturn_ref.*attenuation;
            trfsin=trfunc.*sin(current_targetangle);
            trfcos=trfunc.*cos(current_targetangle);
            
%             jjj=isnan(theta);
%             if jjj>0; disp([ 'help  '  num2str(jj) '    ' num2str(receiver_range) '    '  num2str(theta) ]);end;
            transfer_function(jj,:)=transfer_function(jj,:) +  trfunc;
            transfer_function_sin(jj,:)=transfer_function_sin(jj,:) +  trfcos;
            transfer_function_cos(jj,:)=transfer_function_cos(jj,:) +  trfsin;
            
        end;
    end;
end;% jj
Amp_tot=abs(transfer_function);
Amp_sin=abs(transfer_function_sin);
Amp_cos=abs(transfer_function_cos);
TL_tot= -20*log10(abs(transfer_function));
TL_sin=-20*log10(abs(transfer_function_sin));
TL_cos=-20*log10(abs(transfer_function_cos));

%Convert to time domain and calculate the pulse respones for  ranges  given by  para.range_phone
N_phones=length(para.range_phone);
time=(1:N)/fs;
for jj=1:N_phones
    receiver_range=para.range_phone(jj);
    signal_spectrum=source_spectrum.*transfer_function(jj,:);
    transfer_spectrum_sin=source_spectrum.*transfer_function_sin(jj,:);
    transfer_spectrum_cos=source_spectrum.*transfer_function_cos(jj,:);
    signal(jj,:)=2*real(ifft(signal_spectrum,N));
    signal_sin(jj,:)=2*real(ifft(transfer_spectrum_sin,N))/rhoc;
    signal_cos(jj,:)=2*real(ifft(transfer_spectrum_cos,N))/rhoc;
    t_red(jj,:)=t;
    % t_real(jj,:)=t+(receiver_range-para.range_source)/c_red;
    
end; %jj over receiver positions
%% plotting the time responses

figure(200);clf;
lines=1.5; fonts=16;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
scale=500;

for jj=1:1: N_phones
    hold on;
    plot(t_real(jj,:),signal(jj,:)*scale+bias(jj)/1000,'k','linewidth',lines)
end;
grid;
xlabel('Time – s'); ylabel ('Range – km');
x_min=0; x_max=5; y_min=-1.0; y_max=5;
axis([x_min x_max y_min y_max]);box on;
title([env.title ':  Sound pressure' ]);
text( xmin, ymax-.1, [ ' Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])



figure(201);clf;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
scale=scale*rhoc;

for jj=1:1: N_phones
    hold on;
    plot(t_real(jj,:),signal_sin(jj,:)*scale+bias(jj)/1000,'r','linewidth',lines);
end;
grid;
xlabel('Time – s'); ylabel ('Range – km');
axis([x_min x_max y_min y_max]);box on
title([env.title ':  Vertical particle velocity' ]);
text( xmin, ymax-.2, [ ' Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])



figure(202);clf;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
for jj=1:1: N_phones
    hold on;
    plot(t_real(jj,:),signal_cos(jj,:)*scale+bias(jj)/1000,'g', 'linewidth',lines) ,
end;
grid;
xlabel('Time – s'); ylabel ('Range – km');
axis([x_min x_max y_min y_max]);box on
title([env.title ':  Horizontal particle velocity' ]);
text( xmin, ymax-.1, [ ' Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m '])


break

figure(2); clf; fonts=16; lines=1;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
scale=.1;
for m=1:1:N_phones
    receiver_range=para.range_phone(m);
    plot(scale*signal(m,:)+ receiver_range/1000,'r', 'linewidth',lines);hold on;
end
xmin=0; xmax=400; ymin=0; ymax=1.2*max(para.range_phone)/1000;
axis([xmin xmax ymin ymax]);grid;
title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m',  ]);
legend('Total pressure')


figure(3); clf; fonts=16; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
scale=20;
for m=1:1:N_phones
    receiver_range=phone_positions(m);
    plot(time,scale*signal_sin(m,:)+ receiver_range/1000,'b');hold on;
end
xmin=0; xmax=max(time)/4; ymin=0; ymax=max(phone_positions)/1000;


figure(4); clf %

fonts=16; lines=1.5;
set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);

scale=20;
for m=1:1:N_phones
    receiver_range=phone_positions(m);
    plot(t,scale*signal_sin(m,:)+ receiver_range/1000,'g');hold on;
end
xmin=0; xmax=max(t)/4; ymin=0; ymax=max(phone_positions)/1000;

axis([xmin xmax ymin ymax]);grid
%plot_time_signals;



% figure(2); clf %
%
% fonts=16; lines=1.5;
% set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
%
% % plotting time signal
% range_no= 50:50:200;
% scale=200;
% for m=20:20:100;
%     current_range=range(m);
%     sig= real(ifft(source_spectrum.*source_spectrum, N/2));
%     vert_signal(m,:)=2*real(ifft(source_spectrum.*Amp_sin(m,:), N));
%     plot(time, scale*vert_signal(m,:)+ current_range);hold on
% end

% figure(3)
% plot(time,source_signal,'r',time,real(sig), 'b')
% axis([ 0 max(time)/5, -4 4])
% legend('Source', 'signal')
%
% break
% x_min=min(range); x_max=para.Rmax/1000;
% dB_min=20; dB_max=120;
% plot(range,Amp_tot,'--r',range, Amp_sin,'b', range,Amp_cos,'g','linewidth',lines);  grid on
%
% %axis([x_min x_max dB_min dB_max])
% legend('Inline', 'Vertical', 'Horizontal');
% xlabel('Range – km');
% ylabel('Amplitude ');
% title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m, ','Freq. = ' num2str(frequency) ' Hz']);

% figure(3); clf %
% %disp('Figure 10 plots tranmission loss on a logaritmic  range scale')
%
% set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
% fonts=16; lines=1.5;
% set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
% x_min=-90; x_max=90; y_min=-90; y_max=90;
%
% plot(saved_angle1,saved_angle2, saved_angle1,saved_angle1, '--r','linewidth',lines);  grid on
%
% % axis([x_min x_max y_min y_max])
% legend('Target angle','Start angle' );
% xlabel('Start angle');    ylabel('Target angle ');
% title([ env.title  ':  Sd=' num2str(z_source) ' m, ','Rd=' num2str(z_receiver) ' m, ']);


