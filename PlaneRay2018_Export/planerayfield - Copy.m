function eigen = planerayfield (HIST, COUNT, env, para)
% planerayfield with bug fix 2017 by Tor Arne Reinen in function selectarrival: Marked TAR 2017 on 2 lines

%planerayfield   Computes the the ray amplitude.

%   Authors: Jens M. Hovem and Shefeng Yan
%   Copyright 2008 Acoustic Research Center, NTNU
%   $Revision: 2.0 $  $Date: 2008/04/12 16:13:25 $

phone_positions = para.range_receiver;
z=env.ssp(:,1);
c=env.ssp(:,2);
c0=interp1(z,c,para.z_source,'linear');
c_target=interp1(z,c,para.z_receiver,'linear');

K2 = ceil(size(HIST.THETA,2)/2);
K2 = min(K2,para.countlength);

class=HIST.class;
range_int=phone_positions;
Nr = length(range_int);
NKlass=max(max(class));
ZZ=ones(Nr, NKlass*K2);
rayamp=ZZ*0; distance=ZZ*inf; delay=nan(size(ZZ)); startangle=nan(size(ZZ)); count=nan(7,NKlass*K2);
cnt=0;
for Klass=1:NKlass
    if Klass>1 && Klass<=17
        K3=K2;
    else
        K3=1;
    end;        
    for ii=1:K3
      
        if Klass<=17
            categ = category (Klass,ii);
        else
            categ_acs = [COUNT.BOTTOM(class==Klass),COUNT.SURFACE(class==Klass),COUNT.TURNING_ABOVE(class==Klass),COUNT.TURNING_BELOW(class==Klass),sign(HIST.THETA(class==Klass))];
            categ = categ_acs(1,:);
        end;
        [T_temp, R_temp, DIS_temp, A_temp, TA_temp] =selectarrival( categ(1),categ(2),categ(3),categ(4),categ(5), HIST,COUNT);
        T_temp=T_temp(R_temp>0); A_temp=A_temp(R_temp>0); TA_temp=TA_temp(R_temp>0); R_temp=R_temp(R_temp>0); DIS_temp=DIS_temp(R_temp>0);
        if ~isempty(A_temp)
            cntold=cnt; 
            if issorted(R_temp) || issorted(-R_temp)
                [amp_int,dis_int,A_int, TA_int,T_int]=getamp(phone_positions,R_temp,DIS_temp,A_temp,TA_temp,T_temp,c0, c_target);
                cnt=cnt+1;
            else % double eigenray
                clear amp_int;clear dis_int;clear A_int;clear TA_int;clear T_int;
                peakind=find(diff(sign(diff(R_temp))))+1; peakind=[1,peakind.',length(R_temp)];
                for ind=1:length(peakind)-1
                    [amp_int(ind,:),dis_int(ind,:),A_int(ind,:), TA_int(ind,:),T_int(ind,:)]=getamp(phone_positions,R_temp(peakind(ind):peakind(ind+1)),DIS_temp(peakind(ind):peakind(ind+1)),...
                        A_temp(peakind(ind):peakind(ind+1)),TA_temp(peakind(ind):peakind(ind+1)),T_temp(peakind(ind):peakind(ind+1)) ,c0, c_target);
                end;
                cnt=cnt+length(peakind)-1;
            end;
            rayamp(:,cntold+1:cnt)     = amp_int.';
            distance(:,cntold+1:cnt)   = dis_int.';
            startangle(:,cntold+1:cnt) = A_int.';
            targetangle(:,cntold+1:cnt) = TA_int.';
          
            delay(:,cntold+1:cnt)      = T_int.';  
            count(:,cntold+1:cnt)      = [categ,Klass,ii].'*ones(1,cnt-cntold);
        end;
    end %ii
end %Klass

rayamp      = rayamp(:,1:cnt);
distance    = distance(:,1:cnt);
startangle  = startangle(:,1:cnt);
targetangle = targetangle(:,1:cnt);
delay       = delay(:,1:cnt);
count       = count(:,1:cnt);

eigen.startangle=startangle;
eigen.targetangle=targetangle;
eigen.timedelay=delay;
eigen.rayamp=rayamp;
eigen.distance=distance;
eigen.count=count;
eigen.range=phone_positions;
save -v7.3 eigen eigen


function[amp_int,distance_int,A_int, TA_int,T_int]=getamp(range_int,R,DIS,A,TA,T ,c0, c_target)
if all(abs(diff(A))<5*eps)
    mid=round((length(A)+1)/2);
    A  = A(mid);
    TA = TA(mid);
end;
if length(A)>=2
    method='linear';
    last1= nan;
    T_int=interp1(R,T,range_int, method,last1);
    A_int=interp1(R,A,range_int, method,last1);
    TA_int=interp1(R,TA,range_int, method,last1);
    distance_int=interp1(R,DIS,range_int, method,last1);
    [amp TL]=ray_loss(R, A*pi/180, TA*pi/180,c0, c_target);
    amp_int=interp1(R(1:length(TL)), amp,range_int, method,0);
elseif length(A)==1
    if abs(A)<10*eps
        T_int    = range_int/c_target;
        A_int    = A*ones(size(range_int));
        TA_int   = TA*ones(size(range_int));
        amp_int  = 1./range_int;
        distance_int=range_int;
    else
        T_int    = nan;
        A_int    = nan;
        TA_int   = nan;
        amp_int  = 0;
        distance_int=nan;
    end;
end;


function[amp,TLdB]=ray_loss (range,theta_0, theta_target,c0,c_target);
% Computes the the transmission loss in dB (TLdB) and the ray amplitude (amp)
% for distances given by range, initial angle (theta_0), angle at receiver theta_target).
% c0 and c_target are the sound speed at the source and the receiver depths respectively

% range in meter
% angles in radians
%c0 = sound speed at source depth
%c_target = sound speed at receiver depth
%theta_0 = initial angle in radians
%theta_target = angle at the target in radians
% filterlength=length of smoothing filter;

filterlength=1;
X=1/filterlength*ones(1,filterlength);
R=range;
A=theta_0; 
TA=theta_target;
w=diff(A)./diff(R);
if length(w)==1
    neww=[w;w];
else
    newR=mean([R(1:end-1),R(2:end)],2);
    neww=interp1(newR,w,R,'linear','extrap');
    neww([1,end])=neww([1,end]).*[mean(sign(neww([1,2])));mean(sign(neww([end-1,end])))];
end;
I=(1./R).*(cos(A)./(sin(TA)+eps)).*neww;
% I=filter(X,1,I);
I=abs(I);
amp=sqrt(I);amp(amp>1)=eps;
TLdB=-10 *log10(abs(I)+10^(-20));
TLdB(TLdB<0 | TLdB>=200)=inf;


function[time,range,distance,angle,angle_target,label]=selectarrival(bottom,surface,turning_above,turning_below,direction,HIST,COUNT)

eig_range=HIST.RANGE;
eig_distance=HIST.DISTANCE;
eig_theta=HIST.THETA;
eig_angle=HIST.ANGLE;
eig_time =HIST.TIME;
surface_count=COUNT.SURFACE;
bottom_count=COUNT.BOTTOM;
turning_above_count=COUNT.TURNING_ABOVE;
turning_below_count=COUNT.TURNING_BELOW;
eig_theta=eig_theta*180/pi;
eig_angle=eig_angle*180/pi;

pick=[bottom, surface, turning_above, turning_below];
label=(bottom_count==pick(1) & surface_count==pick(2) & turning_above_count==pick(3) & turning_below_count==pick(4));
time=eig_time(label);time=time(~isnan(time));
range=eig_range(label);range=range(~isnan(range));
distance=eig_distance(label);distance=distance(~isnan(distance));
angle=eig_theta(label);angle=angle(~isnan(angle));
angle_target=eig_angle(label);angle_target=angle_target(~isnan(angle_target));
if direction>0;
    distance=distance(angle>0);  % TAR 2017 
    time=time(angle>0); range=range(angle>0); angle_target=angle_target(angle>0); angle=angle(angle>0);
    label(label>0 & eig_theta<0)=0;

end;
if direction<0;
    distance=distance(angle<0);  % TAR 2017 
    time=time(angle<0); range=range(angle<0); angle_target=angle_target(angle<0); angle=angle(angle<0);
    label(label>0 & eig_theta>0)=0;    
end;
