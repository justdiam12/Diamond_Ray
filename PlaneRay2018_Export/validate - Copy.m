function [ eigresult] = validate( para, env, SUFBOT,eigen, threshold, plots)
%%%%  Date: 2015/03/25
%% Modified by intrducing a threshold value for  eigenrays

%% This program analyses the properties of the eigenrays to one or more specific point.
%% The ranges to these points are given in the parameter eigrange in the input program.
%% If eigrange is not set you will be asked to input the range
%% Plots only eigen rays with aplitude greater than the value of threshold,
%% but if threshold is not set the default value is zero

%% The output of the program is contained in the struct eigresult
%% eigresult.range
%% eigresult.angle
%% eigresult.amplitude
%% eigresult.time
%% eigresult.targetangle
%% eigresult.bottom_range]
%% eigresult.bottom_angle
%% eigresult.surface_range
%% eigresult.surface_angle
if nargin< 6; plots=0;end;
if nargin< 5; threshold=0;end;

cp1=env.cp1; cp2=env.cp2;cs2=env.cs2;rho1=env.rho1;rho2=env.rho2;d=env.laythick;
frequency=para.carrier_frequency+ para.frequency(1);
lambda=env.cp1/frequency;
f=para.frequency(1); c_depth=1500;

eigpara=para;
rayamp=eigen.rayamp;
eigenangle=eigen.startangle;
range=eigen.range;%
c_red=para.c_red;
timedelay=eigen.timedelay;
target_angle=eigen.targetangle;
phone_positions=para.range_receiver;
phone_space=mean(phone_positions(2:end)-phone_positions(1:end-1));
if isempty(para.eigrange)
    receiver_range_str=input(['Range for eigenray calculations [ ',num2str(min(phone_positions)),....
        ' : ',num2str(phone_space),' : ',num2str(max(phone_positions)),' ] m =  '],'s');
    receiver_range_input = str2num(receiver_range_str);
else
    receiver_range_input=para.eigrange;
end;
desired_ranges= receiver_range_input;
N_ranges=length(desired_ranges);
for n=1:N_ranges
    select_range= desired_ranges(n);
    if select_range> para.Rmax;  disp('Receiver range too large'); break;   end;
    [temp, jj] = min(abs(range-select_range));
    select_ang=eigenangle(jj,:);
    receiver_angle=target_angle(jj,:);
    time=timedelay(jj,:);
    amp=rayamp(jj,:);
    k=~isnan(time);
    time=time(k);
    amp=amp(k);
    receiver_angle=receiver_angle(k);
    select_ang=select_ang(k);
    
    
    time_red=time-select_range/para.c_red;
    ind_selectang= ~isnan(select_ang);
    select_ang=select_ang(ind_selectang);
    
    select_angle=select_ang;
    eigpara.scan=select_angle;
    select_angle_r=fix(100*select_angle)/100;
    amplitude=amp(ind_selectang);
    amplitude=amplitude(amplitude~=0);
    
    if ~para.quiet
        for i=1:min([length(select_ang) length(time_red)])
            disp(['Eigenangle =' num2str(select_angle_r(i)) '(deg)    '  'Amp-rel = ' num2str(amplitude(i)) ,......
                ' Travel time = ' num2str(1000*time_red(i)) ' ms'])
        end;
    end;
    
    disp(' Plot the eigenray traces and the envibonment ')
    eigpara.start_angles=select_angle;
    figure(plots)
    [ HIST_temp COUNT_temp SUFBOT_temp Rays] = tracerays(env, eigpara, select_range);
     hold on; plot(select_range/1000,para.z_receiver,'ro');
        % End plotting eigenrays
 
    % Continue the calculations of eigenray properies
    save  SUFBOT_temp SUFBOT_temp
    bottom_angle=SUFBOT_temp.BOTTOM_ANGLE;
    bottom_range=SUFBOT_temp.BOTTOM_RANGE;
    surface_angle=SUFBOT_temp.SURFACE_ANGLE;
    surface_range=SUFBOT_temp.SURFACE_RANGE;
    
    
    k=~isnan(time);  time=time(k); N=length(time);
    eigresult.range=select_range;
    eigresult.angle=select_angle;
    eigresult.amplitude=amplitude;
    eigresult.time=time;
    eigresult.targetangle=receiver_angle(k);
    eigresult.bottom_range=bottom_range;
    eigresult.bottom_angle=bottom_angle;
    eigresult.surface_range=surface_range;
    eigresult.surface_angle=surface_angle;
    
    select_range=eigresult.range;
    select_angle=eigresult.angle;
    
    eigpara=para;
    rayamp=eigen.rayamp;
    eigenangle=eigen.startangle;
    range=eigen.range;
    Ni=length(select_angle);
    bottom_ref=ones(1,Ni);sgn=ones(1,Ni);
    for ii=1:Ni;
        sur_range=surface_range(ii,:);
        sur_range=sur_range(sur_range<select_range);
        sgn(ii)=(-1)^(length(sur_range));
        range_bottom= bottom_range(ii,:);
        theta_bottom =bottom_angle(ii,:);
        range_bottom=range_bottom(range_bottom<select_range);
        theta_bottom=theta_bottom(range_bottom<select_range);
        
        Nk=length(theta_bottom);
        ref=1;
        for kk=1:Nk;
            angle=theta_bottom(kk);
            R_rough= R_coh( angle,env.sigma_bottom, lambda );
            temp= R_rough.*abs(R_bottom(angle,f,c_depth,cp1,cp2,cs2,rho1,rho2,d));
            ref=ref*temp;
        end;
        bottom_ref(ii)=ref;
    end
    new_amp=amplitude.*sgn.*bottom_ref;
    amplitude=new_amp;
    eigresult.amplitude=amplitude;
    
    
    if isempty(select_angle); disp( 'No eigenangle to the selected range')
    else
        % Plot only eigen rays with amplitude greater than the value of threshold
        theta_select=eigresult.angle(abs(eigresult.amplitude)>threshold);
        eigpara.start_angles=theta_select;
        save eigpara eigpara;
        save  SUFBOT_temp SUFBOT_temp ;   save eigresult   eigresult;
        plots=1;
        if plots>0;
            
            figure(100+(n-1)*5); clf;  fonts=18; lines=3;
            bw=.3;
            set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
            plot(select_angle,new_amp,'bo','linewidth',lines);
            if N>1; hold on;   bar(select_angle,new_amp,bw);end;
            xmin=-90;xmax=90; ymax=max(new_amp)*1.1; ymin=-ymax;
            axis( [xmin xmax ymin ymax]);
            ylabel('Ray amplitue ');xlabel('Eigenangles – deg') ;grid
            
            title([env.title ':  Sd:' num2str(para.z_source) ' m  Rd: ' ...
                num2str(para.z_receiver) ' m Range: ', num2str(select_range/1000) ' km']);
           
            figure(100+(n-1)*5+1);clf;
            lines =3;
            set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
            
            time_min=min(time_red);
            time_ms=(time_red-time_min)*500;
            N=length(time_ms);
            
            x=time_ms; y=new_amp;
            plot(x, y(1:N),'ro', 'linewidth', lines);
            OK=1;
            %Don't use bar if duplicated x values
            %if any(x(2:end)-x(1:end-1) == 0); OK=0; end;
         %   if OK==1; hold on; bar(time_ms,   new_amp(1:N));end;
            xlabel( 'Reduced time delay – ms')
            ylabel('Ray amplitue ');grid
            xmax=0.5*max(time_ms);xmin=-xmax;
            xmin=-0.1*max(time_ms); 
         axis( [xmin xmax ymin ymax]);
            
         %   title([env.title ':  Sd:' num2str(para.z_source) ' m  Rd: ' ....
             %   num2str(para.z_receiver) ' m Range: ', num2str(select_range/1000) ' km']);
            
            
            
            
            figure(100+(n-1)*5+2);clf;
            set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
            %bw=0.3;
            targetangle=eigresult.targetangle;
            time_min=min(time_red);
            time_ms=(time_red-time_min)*1000;
            N=length(time_ms);
            plot(targetangle,   new_amp(1:N),'ro', 'linewidth', lines);
            if N>1; hold on; bar(targetangle,   new_amp(1:N),bw);end;
            xlabel( 'Target angle  – deg')
            ylabel('Ray amplitue ');grid
           xmin=-90;xmax=90; ymax=max(new_amp)*1.1; ymin=-ymax;
            axis( [xmin xmax ymin ymax]);
           % title([env.title ':  Sd:' num2str(para.z_source) ' m  Rd: '  num2str(para.z_receiver) ' m Range: ', num2str(select_range/1000) ' km']);
            figure(100+(n-1)*5+4); fonts=16;lines=2;
            set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
            
            figure(100+(n-1)*5+3);clf;
            set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
            bw=1;
            time_min=min(time_red);
            time_ms=(time_red-time_min)*1000;
            N=length(time_ms);
            semilogy(select_angle,time_ms,'bo','linewidth',lines);grid
            
            if N>2; hold on; bar(select_angle,   time_ms,bw);end;
            ylabel( 'Reduced time delay – ms')
            xlabel('Eigenangles – deg')
           axis([ -90 90 1 1.1*max(time_ms)])
           % title([env.title ':  Sd:' num2str(para.z_source) ' m  Rd: '  num2str(para.z_receiver) ' m Range: ', num2str(select_range/1000) ' km']);
            
        end 
    end
end;


