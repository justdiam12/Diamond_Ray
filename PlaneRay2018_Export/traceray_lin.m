function [R, Z, T, eig, count, sufbot] = traceray_lin(theta_deg, env, para,R_stop)
%  Traces the ray corresponding to a particular take-off angle when the
%   sound speed gradient is zero.
%   Authors: Jens M. Hovem and Shefeng Yan
%   Copyright 2008 Acoustic Research Center, NTNU
%   Revision:  2014/06/15

if nargin==3;R_stop=para.Rmax; end;
ssp          = env.ssp;
bathy        = env.bathy;
bathslope    = env.bathslope;
range_source = para.range_source;
z_source     = para.z_source;
z_receiver   = para.z_receiver;
range_receiver=para.range_receiver;
N            = para.N;
countlength  = para.countlength;
R_max        = para.Rmax;
bottom_stop  = para.bottom_stop;
if size(ssp,1)<size(ssp,2),ssp=ssp.';end;
z     = ssp(:,1);
c     = ssp(:,2);
del_z = z(2)-z(1);
R_bathy   = bathy(1,:);
Z_bathy   = bathy(2,:);
R_slope   = bathslope(1,:);
ang_slope = bathslope(2,:);
countbuff = 0;

% Set up start values
target_count=0; bottom_hits=0; surface_hits=0; turning_above=0; turning_below=0;
eig_theta=zeros(1,countlength); eig_time=zeros(1,countlength); eig_range=zeros(1,countlength);eig_angle=zeros(1,countlength);
eig_distance=zeros(1,countlength);
surface_count=zeros(1,countlength); bottom_count=zeros(1,countlength); turning_above_count=zeros(1,countlength); turning_below_count=zeros(1,countlength);
surface_range=zeros(1,countlength); bottom_range=zeros(1,countlength); bottom_angle=zeros(1,countlength);
surface_angle=zeros(1,countlength);
delta_t=0; delta_x=0; delta_z=0;
k_target=round(z_receiver/del_z)+1; k_start=round(z_source/del_z)+1;k=k_start;

% Initial settings for the ray starting with theta_deg
theta_0=theta_deg*pi/180; theta_current=theta_0;
theta_in=theta_0;
dir=sign(theta_deg);
x_temp=range_source; z_temp=z_source; t_temp=0;
r_temp=0;
k_temp=k_target;
n_temp=0;

if abs(theta_deg)<5*eps % theta_deg = 0
    nreceiver = length(range_receiver);
    R = [0,range_receiver,R_max];
    Z = z_source*ones(1,nreceiver+2);
    T = R/c(k_start);
    if k_target==k_start % target depth has been detected when theta_deg = 0
        eig_range    = range_receiver([1,min(2,length(range_receiver)),end]);
        eig_distance = eig_range;
        target_count = length(eig_range);
        eig_time     = eig_range/c(k_start);
        eig_theta    = theta_0*ones(1,target_count);
        eig_angle    = eig_theta;
        bottom_count = zeros(1,target_count);
        surface_count= zeros(1,target_count);
        turning_above_count = zeros(1,target_count);
        turning_below_count = zeros(1,target_count);
        surface_hits = 0;
        bottom_hits  = 0;
    end;
else
    R=zeros(1,N);T=zeros(1,N);Z=zeros(1,N);
    ibathy=1; islope=1;
    for n=1: N %Counter along the ray
        if x_temp>=R_stop || x_temp<0 ||  abs(theta_current)>=pi/2, break; end;
        p_current=cos(theta_in)./c(z==z_source);
        if x_temp>=R_max || x_temp<0 || abs(theta_current)>=pi/2, break; end;
        while x_temp>=R_bathy(ibathy+1)
            ibathy=ibathy+1;
        end;
        depth_current=(Z_bathy(ibathy)+Z_bathy(ibathy+1))/2;
        
        while x_temp>=R_slope(islope+1)
            islope=islope+1;
        end;
        angle_current=ang_slope(islope);
        if dir==1; % Down going wave
            k_bottom=fix(depth_current/del_z);
            if k>k_bottom; %Bottom detected;
                if bottom_stop>0;
                    if bottom_hits>=bottom_stop-1; break; end;
                end
                
                if bottom_stop~=0;
                    if bottom_hits>=bottom_stop-1; break; end;
                end
                bottom_hits=bottom_hits+1; dir=-1; %Increment number of bottom hits and turn
                delta_x=0; delta_z=0;  delta_t=0;
                k_bottom=fix(depth_current/del_z); c_bottom=c(k_bottom);
                
                theta_in=acos(p_current*c_bottom); theta_out=theta_in+2*angle_current;
                theta_in=theta_out;
                theta_current=acos(p_current*c_bottom)+2*angle_current;
                
                
                if imag(theta_current)~=0,break;end;
                if theta_current>=pi/2,break;end;
                p_current=cos(theta_out)/c_bottom;
                bottom_range(bottom_hits)=x_temp;
                bottom_angle(bottom_hits)=theta_current-angle_current;
                
                if theta_current<0 % continue down going
                    dir=1;
                    delta_z=delta_z+del_z;
                    delta_x=delta_x+del_z/abs(tan(theta_current));
                    delta_t=delta_t+del_z/abs((sin(theta_current))*c_bottom);
                    k=k+1;
                else
                    newibathy=ibathy;
                    depth_current=(Z_bathy(newibathy)+Z_bathy(newibathy+1))/2;
                    while z_temp+delta_z>=depth_current % make sure the ray return to water
                        delta_z=delta_z-del_z;
                        delta_x=delta_x+del_z/abs(tan(theta_current));
                        delta_t=delta_t+del_z/abs((sin(theta_current))*c(k));
                        k=k-1;
                        
                        if newibathy+1>=length(R_bathy),break;end;
                        depth_current=(Z_bathy(newibathy)+Z_bathy(newibathy+1))/2;
                    end;
                end;
            else
                if k<length(c)
                    
                    
                    if  p_current*c(k+1)<1 % No turning in next layer below
                        ck=c(k); ck1=c(k+1);
                        theta_current=acos(p_current*ck1);
                        delta_x=del_z/abs(tan(theta_current));
                        delta_t=del_z/abs((sin(theta_current))*ck1);
                        delta_z=del_z;
                        k=k+1;
                    else % Turning point within the next layer below
                        turning_below=turning_below+1; dir=-dir;
                        delta_x=0;delta_t=0;delta_z=0;
                    end;% End turning point
                end; % End down going wave
            end
        else  % Up going wave
            
            if k==1 % surface hit
                dir =1;
                surface_hits=surface_hits+1;
                surface_range(surface_hits)=x_temp;
                theta_current=acos(p_current*c(1));
                surface_angle(surface_hits)=theta_current;
                delta_x=0;
                delta_z=0;  delta_t=0;
                k=1;
            elseif z_temp>=depth_current; %Bottom detected
                if bottom_stop==1 ; break; end;
                if n-n_temp==1;break;end;n_temp=n;
                bottom_hits=bottom_hits+1; %Increment number of bottom hits
                while x_temp-real(delta_x)<R_bathy(ibathy) && round((z_temp-delta_z)/del_z)-round(Z_bathy(ibathy)/del_z)>=0
                    ibathy=ibathy-1;
                end;
                if p_current*c(k+1)>=1 || p_current*c(k)>=1
                    theta_current=0;
                    zm=0;
                    delta_z=0;
                    delta_x=R_bathy(ibathy)-x_temp;
                    while x_temp+delta_x<R_slope(islope)
                        islope=islope-1;
                    end;
                    angle_current=ang_slope(islope);
                else
                    
                    if imag(theta_current)~=0,break;end;
                    needdo=1;
                    while needdo
                        needdo=0;
                        angle_current=ang_slope(islope);
                        zm=(R_bathy(ibathy+1)-(x_temp-real(delta_x)))/(1/tan(theta_current)-1/tan(angle_current));
                        zm=min(zm,del_z);
                        delta_x_tmp=-delta_x+zm/tan(theta_current);
                        while x_temp+delta_x_tmp<R_slope(islope)
                            islope=islope-1;
                            needdo=1;
                        end;
                    end;
                    delta_x=-delta_x+zm/tan(theta_current);
                    delta_z=-delta_z-zm;
                end;
                
                k_bottom=fix(depth_current/del_z);
                delta_t=-abs(delta_x/abs(cos(theta_current)*c(k_bottom)));
                
                x_temp=x_temp+real(delta_x);
                z_temp=z_temp+delta_z;
                r_temp=r_temp-sqrt(delta_x^2+delta_z^2);
                t_temp=t_temp+delta_t;
                T(n-1)=real(t_temp);
                Z(n-1)=z_temp;
                R(n-1)=x_temp;
                c_bottom=c(k);
                theta_current=dir*acos(p_current*c_bottom)+2*angle_current;
                if imag(theta_current)~=0,break;end;
                
                p_current=cos(theta_current)/c_bottom;
                bottom_range(bottom_hits)=x_temp;
                bottom_angle(bottom_hits)=theta_current-angle_current;
                theta_current=theta_current*dir;
                
                delta_z=-del_z+zm;
                delta_x=-delta_z/tan(abs(theta_current));
                delta_t=abs(delta_x/abs(cos(theta_current)*c(k)));
                delta_x=delta_x+del_z/abs(tan(theta_current));
                delta_z=delta_z-del_z;
                delta_t=delta_t+del_z/abs(sin(theta_current)*c_bottom);
                k=k-1;
            else  %No surface hit this time
                ck=c(k); ck1=c(k+1);
                k=k-1;
                g=(c(k+1)-c(k))/(z(k+1)-z(k));
                
                if p_current*c(k)<1 % No turning within next layer above
                    theta_current=acos(p_current*ck1);
                    delta_x=del_z/tan(theta_current);
                    delta_t=del_z/(sin(theta_current)*ck1);
                    delta_z=-del_z;
                else
                    % Turning point within next layer above
                    turning_above=turning_above+1; dir=-dir;
                    delta_z=0; delta_x=0; delta_t=0; k= k+1;
                end; % End turning point
            end;% Detection of surface hits
        end % Updating ray position and traveltime.
        
        
        x_temp=x_temp+real(delta_x);
        z_temp=z_temp+delta_z;
        r_temp=r_temp+sqrt(delta_x^2+delta_z^2);
        t_temp=t_temp+delta_t;
        % End updating ray position and traveltime.
        
        % Determination of range to target(receiver) depth.
        if k==k_target
            target_countold=target_count;
            target_count=target_count+1; % target_count is the number of times target depth has been detected,
            eig_range(target_count)=x_temp;
            eig_distance(target_count)=r_temp;
            eig_time(target_count)=t_temp;
            eig_theta(target_count)=theta_0;
            eig_angle(target_count)=dir*theta_current;
            bottom_count(target_count)=bottom_hits;
            turning_above_count(target_count)=turning_above;
            turning_below_count(target_count)=turning_below;
            surface_count(target_count)=surface_hits;
            countbuff=1;
        end;
        if countbuff
            if k_target==k_temp %when the ray trace is horizontal.
                target_count = target_countold;
            else
                target_count = target_countold+1;
                countbuff=0;
            end;
        end;
        
        T(n)=real(t_temp);
        Z(n)=z_temp;
        R(n)=x_temp;
        k_temp=k;
    end % n loop
    R=R(R>0);
    Z=Z(R>0);
    R = [0,       R];
    Z = [z_source,Z];
    T = [0,       T];
end;
count.targethit     = target_count;
count.surfacehit    = surface_hits;
count.bottomhit     = bottom_hits;
count.bottom        = bottom_count(1:target_count);
count.surface       = surface_count(1:target_count);
count.turning_above = turning_above_count(1:target_count);
count.turning_below = turning_below_count(1:target_count);

eig.angle = eig_angle(1:target_count);
eig.range = eig_range(1:target_count);
eig.distance=eig_distance(1:target_count);
eig.time  = eig_time(1:target_count);
eig.theta = eig_theta(1:target_count);

bottom_angle=bottom_angle(1:bottom_hits);
bottom_range=bottom_range(1:bottom_hits);
surface_angle=surface_angle(1:surface_hits);
surface_range=surface_range(1:surface_hits);
sufbot.bottom_angle=bottom_angle;
sufbot.bottom_range=bottom_range;
sufbot.surface_angle=surface_angle;
sufbot.surface_range=surface_range;

if ~isreal(eig.time)
    warning([mfilename,' doesn''t work at this ssp']);
end;