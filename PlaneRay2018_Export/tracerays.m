function [HIST, COUNT, SUFBOT, Rays] = tracerays(env, para, R_stop)
%   Traces the ray corresponding to a particular start angle 
%   Authors: Jens M. Hovem and Shefeng Yan
%   Date 2015/05/15

if nargin==2; R_stop=para.Rmax; end;

if ~para.quiet,disp('Start ray tracing ...');tic;end
plotenv(env,para);
start_theta = para.start_angles;

nrays=length(start_theta);
count_targethit=zeros(1,nrays);
count_surfacehit=zeros(1,nrays);
count_bottomhit=zeros(1,nrays);
count_raylength=zeros(1,nrays);
countlength=para.countlength;

HIST_TIME=nan(nrays,para.countlength);
HIST_RANGE=nan(nrays,para.countlength);
HIST_THETA=nan(nrays,para.countlength);
HIST_ANGLE=nan(nrays,para.countlength);
HIST_DISTANCE=nan(nrays,para.countlength);

BOTTOM_COUNT=zeros(nrays,para.countlength);
SURFACE_COUNT=zeros(nrays,para.countlength);
TURNING_ABOVE_COUNT=zeros(nrays,para.countlength);
TURNING_BELOW_COUNT=zeros(nrays,para.countlength);
% 



SUFBOT.BOTTOM_ANGLE=ones(nrays,para.countlength)*1000;
SUFBOT.BOTTOM_RANGE=2*para.Rmax*ones(nrays,para.countlength);
SUFBOT.SURFACE_ANGLE=nan(nrays,para.countlength);
SUFBOT.SURFACE_RANGE=nan(nrays,para.countlength);
ALLX=nan(nrays,para.N);
ALLZ=nan(nrays,para.N);

for ray_number=1:nrays
    theta_deg=start_theta(ray_number);
    if ~para.quiet
        disp([' Ray no ' num2str(ray_number) '  of  ' num2str(nrays) ' Angle =  ' num2str(theta_deg) '  degs']); end;
    switch para.rayopt
        case 1
            [ X, Z, TT, eig, count, sufbot] = traceray_lin (theta_deg, env, para,R_stop);
        case 0
            [ X, Z, TT, eig, count, sufbot] = traceray_circ(theta_deg, env, para,R_stop);
    end;
    X=X(Z>=0); Z=Z(Z>=0);
   plot(X/1000,Z,'k','linewidth',0.5);
    
    
    count_targethit(ray_number) =min(count.targethit,countlength);
    count_surfacehit(ray_number)=count.surfacehit;
    count_bottomhit(ray_number) =count.bottomhit;
    
    HIST_TIME (ray_number,1:min(length(eig.time), countlength))      =eig.time (1:min(length(eig.time),countlength));
    HIST_RANGE(ray_number,1:min(length(eig.range),countlength))      =eig.range(1:min(length(eig.range),countlength));
    HIST_THETA(ray_number,1:min(length(eig.theta),countlength))      =eig.theta(1:min(length(eig.theta),countlength));
    HIST_ANGLE(ray_number,1:min(length(eig.angle),countlength))      =eig.angle(1:min(length(eig.angle),countlength));
    HIST_DISTANCE(ray_number,1:min(length(eig.distance),countlength))=eig.distance(1:min(length(eig.distance),countlength));
    
    BOTTOM_COUNT(ray_number,1:length(count.bottom))=count.bottom;
    SURFACE_COUNT(ray_number,1:length(count.surface))=count.surface;
    TURNING_ABOVE_COUNT(ray_number,1:length(count.turning_above))=count.turning_above;
    TURNING_BELOW_COUNT(ray_number,1:length(count.turning_below))=count.turning_below;
    
    SUFBOT.BOTTOM_ANGLE(ray_number,1:length(sufbot.bottom_angle))=sufbot.bottom_angle;
    SUFBOT.BOTTOM_RANGE(ray_number,1:length(sufbot.bottom_range))=sufbot.bottom_range;
    SUFBOT.SURFACE_ANGLE(ray_number,1:length(sufbot.surface_angle))=sufbot.surface_angle;
    SUFBOT.SURFACE_RANGE(ray_number,1:length(sufbot.surface_range))=sufbot.surface_range;
    
    count_raylength(ray_number)=length(X);
 
        if count_raylength(ray_number) > para.N
            ALLX=[ALLX,nan(nrays,count_raylength(ray_number)-para.N)];
            ALLZ=[ALLZ,nan(nrays,count_raylength(ray_number)-para.N)];
            para.N=count_raylength(ray_number);
        end;
        ALLX(ray_number,1:count_raylength(ray_number))=X;
        ALLZ(ray_number,1:count_raylength(ray_number))=Z;

end;
maxtargethit=max(count_targethit);
maxsurfacehit=max(count_surfacehit);
maxbottomhit=max(count_bottomhit);

HIST.TIME=HIST_TIME(:,1:maxtargethit);
HIST.RANGE=HIST_RANGE(:,1:maxtargethit);
HIST.THETA=HIST_THETA(:,1:maxtargethit);
HIST.ANGLE=HIST_ANGLE(:,1:maxtargethit);
HIST.DISTANCE=HIST_DISTANCE(:,1:maxtargethit);

COUNT.BOTTOM=BOTTOM_COUNT(:,1:maxtargethit);
COUNT.SURFACE=SURFACE_COUNT(:,1:maxtargethit);
COUNT.TURNING_ABOVE=TURNING_ABOVE_COUNT(:,1:maxtargethit);
COUNT.TURNING_BELOW=TURNING_BELOW_COUNT(:,1:maxtargethit);

SUFBOT.BOTTOM_ANGLE =SUFBOT.BOTTOM_ANGLE(:,1:maxbottomhit);
SUFBOT.BOTTOM_RANGE =SUFBOT.BOTTOM_RANGE(:,1:maxbottomhit);
SUFBOT.SURFACE_ANGLE=SUFBOT.SURFACE_ANGLE(:,1:maxsurfacehit);
SUFBOT.SURFACE_RANGE=SUFBOT.SURFACE_RANGE(:,1:maxsurfacehit);
save  SUFBOT SUFBOT

     Rays.X=ALLX(:,1:max(count_raylength));
     Rays.Z=ALLZ(:,1:max(count_raylength));

if ~para.quiet,disp(['End ray tracing. Total time: ',num2str(toc),' s']);end;
