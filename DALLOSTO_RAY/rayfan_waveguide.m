%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RAY FAN
%PROGRAMMER: DAVID DALLOSTO
%3/5/2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Range,Depth,time] = rayfan_waveguide(theta,RANGERES,RR,SD,dcz,cw,zw,cb,zb)
%%%%%%%%%%%%%%%%
% INPUTS
% cw = water soundspeed profile
% zw = soundspeed profile depth index
% cb = seabed soundspeed profile (layered)
% zb = depth of seabed profile
% SD/SR = source depth and range
% numrays = number of rays in rayfan
% minangle/maxangle = angles for ray fan
%%%%%%%%%%%%%%%%%
% OUTPUTS
% Depth = Matrix of Depths
% Range = Matrix of Ranges
% time = Matrix of travel times
%%%%%%%%%%%%%%
maxrange = RR+RANGERES;
tangvector = -dcz./cw; %Tangent Vector to ray V = -1/c*(dc/dz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ntheta = 1:length(theta)

%%%%%%%%%%%%%%%%Initialization%%%%%%%%%%%%%%%%%%%
dR = RANGERES;
thetalaunch = -theta(ntheta);
R(1) = 0;
D(1) = SD;
t(1) = 0;
tau(1) = tand(thetalaunch);  % TAU = tan(theta);


ii=find(zw==SD); %starting depth index
c0=cw(ii);       %forTLtest
s(1)=dR/cos(atan(tau(1)));     %forTLtest
phase(1) = 0;

i=2;maxi =  ceil(maxrange/RANGERES)+1;
while i <= maxi %calculations out to 2000 m (make sure include RR inside number)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Insert R-K 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test=abs(zw-D(i-1));
ztemp=find(test==min(test),1,'first');
k1tau=(1+tau(i-1)^2)*tangvector(ztemp);
k1z=tau(i-1);
test=abs(zw-(D(i-1) + k1z*.5*dR));
nuztemp=find(test==min(test),1,'first');
k2tau=(1+(tau(i-1)+k1tau*dR/2)^2)*tangvector(ztemp);
k2z=tau(i-1)+k1tau*dR/2;
test=abs(zw-(D(i-1) + k2z*.5*dR));
ztemp=find(test==min(test),1,'first');
k3tau=(1+(tau(i-1)+k2tau*dR/2)^2)*tangvector(ztemp);
k3z=tau(i-1)+k2tau*dR/2;
test=abs(zw-(D(i-1) + k3z*dR));
ztemp=find(test==min(test),1,'first');
k4tau=(1+(tau(i-1)+k3tau*dR)^2)*tangvector(ztemp);
k4z=tau(i-1)+k3tau*dR;
 	
D(i)=D(i-1)+(dR/6)*(k1z + 2*k2z + 2*k3z + k4z);
R(i)=R(i-1)+dR;
dt = sqrt(1+tau(i-1)^2)/cw(ztemp)*dR;
t(i) = t(i-1)+dt;

tau(i)=tau(i-1)+(dR/6)*(k1tau + 2*k2tau + 2*k3tau(1) + k4tau);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Reflection Logic
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorlevel = 1e-5; 
%%%%%%%%%%%%%%%%%Rounding Error Correction%%%%%%%%%%%%%%%%%%%
        waterdepth = zw(end);
        tolerance = abs(D(i));
 
        if tolerance<errorlevel
            D(i)=0;
        end
        tolerance = abs(waterdepth-D(i));
        if tolerance<errorlevel
            D(i)=waterdepth;
            

        end

        if D(i)<0

            dR = -D(i-1)/tau(i-1);
              i = i-1;          
        elseif D(i)==0
            tau(i) = -tau(i);

        elseif D(i)>waterdepth

            dR = -(D(i-1)-waterdepth)/tau(i-1);
            i = i-1;
        elseif D(i)==waterdepth

            tau(i) = -tau(i);
            
        else
            dR = RANGERES;

        end
        i = i+1;
end


%Must include step to land at RR precisely.  %Section Added 9/22/2014

%Find points that bracket reciever range
Rplus = find(R>RR,1,'first');
Rminus = find(R<RR,1,'last');
%Logic if this falls short of reciever range
if prod(size(Rplus))~=0
    Rplus = Rminus;
    Rminus = Rminus-1;
end
%Compute slope of line
Rang = atan2((D(Rplus)-D(Rminus)),(R(Rplus)-R(Rminus)));
%Extend line from Rminus point to Reciever Range
D(Rminus+1:end) = D(Rminus)+tan(Rang)*(RR-R(Rminus));
R(Rminus+1:end) = RR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting Terms from Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Depth(:,ntheta) = D(1:maxi-1);
Range(:,ntheta) = R(1:maxi-1);
time(:,ntheta) = t(1:maxi-1);
% S(:,ntheta) = s;
% P(:,ntheta) = p;
% Phase(:,ntheta) = phase;
%indexR(ntheta)=find(R>=RR-floor(RANGERES/2) & R<=RR+ceil(RANGERES/2),1,'first');

end

   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% subplot(2,4,[1 5])
% h = plot(cw,zw);
% set(gca, 'ydir','reverse')  %so now positive axis is downward
% title('Sound Speed Profile')
% subplot(2,4,[2:4 6:8])
% plot(0,SD,'o',[0 maxrange],[0 0],'b',[0 maxrange],[waterdepth waterdepth],'k');
% hold on
% plot(real(Range),real(Depth))
% axis([-5,RR+5,-10 waterdepth+zb(end)+10])
% set(gca, 'ydir','reverse')  %so now positive axis is downward
% xlabel('RANGE (m)')
% ylabel('DEPTH (m)')