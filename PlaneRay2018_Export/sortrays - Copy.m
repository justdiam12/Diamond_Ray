function [eigenangle, HIST, COUNT, SUFBOT, Rays, para] = sortrays(HIST, COUNT, SUFBOT, Rays, env, para)
%% sorttrays Sortout the ray classes

%%  Date: 2015/06/25
%% This program is simplification of an earlier program 
%% Input; HIST,COUNT, para,env; most recent PlaneRay results

%% Output is caustic =     caustics=[caustic_range' caustic_angle' caustic_amp'];
%% figure(plots): Ray history and initial angles of the peaks
%% figure(plots+1)Ray amplitudes of peaks as function of range

if nargin<=5 % sortrays(HIST, COUNT, Rays, env, para) used in ploteigstructure and ploteigenray.
    para=env;env=Rays;Rays=SUFBOT;clear SUFBOT;
    SUFBOT.BOTTOM_ANGLE=[]; SUFBOT.BOTTOM_RANGE=[]; SUFBOT.SURFACE_RANGE=[];
end;
if ~para.quiet; disp(['Start sorting and classifying rays']); end;
[HIST, COUNT, SUFBOT, Rays, para, eigangle]= classifyrays (HIST, COUNT, SUFBOT, Rays, env, para);
if nargout,eigenangle=eigangle;end;

HIST_RANGE=HIST.RANGE;

for ray_number=1:size(Rays.X,1)
    X=Rays.X(ray_number,:);
    X=X(~isnan(X));
    Z=Rays.Z(ray_number,:);
    %Z=Z(~isnan(X));
    range=HIST_RANGE(ray_number,:);
    %range=range(~isnan(range));
   % segind = rayseg(X,range);
end   
   