function [HIST, COUNT, SUFBOT, Rays, para_bak, startangle] = classifyrays (HIST, COUNT, SUFBOT, Rays, env, para)

K2 = ceil(size(HIST.THETA,2)/2);
K2 = min(K2,para.countlength);

class =zeros(size(HIST.THETA));class(isnan(HIST.THETA))=nan;
class2=zeros(size(HIST.THETA));class2(isnan(HIST.THETA))=nan;
Nreceiver=length(para.range_receiver);
startangle=nan(Nreceiver,17*K2);
count=nan(8,17*K2);
cnt=0;
for Klass=1:17
    if Klass>1 && Klass<=17
        K3=K2;
    else
        K3=1;
    end;
    for ii=1:K3
        
        categ = category (Klass,ii);
        [T_temp R_temp A_temp TA_temp label] =selectarrival( categ(1),categ(2),categ(3),categ(4),categ(5), HIST,COUNT);
        class(label & class==0)=Klass;class2(class==Klass & class2==0)=ii;
               if ~isempty(A_temp)
            if all(diff(A_temp)<5*eps)
                A_temp  = A_temp(round((length(A_temp)+1)/2));
            end;
            cntold=cnt;
            if length(A_temp)>=2
                 
                if issorted(R_temp) || issorted(-R_temp)
                    A_int=interp1(R_temp,A_temp,para.range_receiver, 'linear',nan);
                    cnt=cnt+1;
                else % double eigenray
%                     %disp('Double eigenray' )
%                     R_caustic=min(R_temp);
%                     A_caustic= A_temp(R_temp==min(R_temp));
%                     disp(['Warning Caustic at range '  num2str(R_caustic) ' meter'])
%  
%                     caustic_position=  [A_temp R_temp];
%                    save  caustic_position caustic_position;
                    clear A_int;
                    peakind=find(diff(sign(diff(R_temp))))+1; peakind=[1,peakind.',length(R_temp)];
                    for ind=1:length(peakind)-1
                        A_int(ind,:)=interp1(R_temp(peakind(ind):peakind(ind+1)),A_temp(peakind(ind):peakind(ind+1)),para.range_receiver, 'linear',nan);
                    end;
                    cnt=cnt+length(peakind)-1;
                end;
            elseif length(A_temp)==1
                if abs(A_temp)<5*eps
                    A_int = A_temp*ones(size(para.range_receiver));
                else
                    A_int = nan;
                end;
                cnt=cnt+1;
            end;
            startangle(:,cntold+1:cnt) = A_int.';
            [iang,ir]=find(label);
            count(:,cntold+1:cnt)      = [categ,Klass,ii,ir(1)].'*ones(1,cnt-cntold);
        end
    end %ii
end %Klass

categ_acs=[COUNT.BOTTOM(class==0),COUNT.SURFACE(class==0),COUNT.TURNING_ABOVE(class==0),COUNT.TURNING_BELOW(class==0),sign(HIST.THETA(class==0))];

categ_acs=sortrows(categ_acs);
if ~isempty(categ_acs)
    categ_acs=categ_acs([1,any(diff(categ_acs,1,1).')]>0,:);
else
    categ_acs=[];
end;
% startangle_acs=nan(size(categ_acs,1),Nreceiver);
para_bak=para; newpara=para; newpara.quiet=1;
for inewklass=1:size(categ_acs,1)
    [T_temp R_temp A_temp TA_temp label] =selectarrival( categ_acs(inewklass,1),categ_acs(inewklass,2),categ_acs(inewklass,3),categ_acs(inewklass,4),categ_acs(inewklass,5), HIST,COUNT);
    class(label & class==0)=17+inewklass;    class2(class==17+inewklass & class2==0)=1;
    cntold=cnt;
    para.searchrays=1;
    if length(A_temp)>=2
        if issorted(R_temp) || issorted(-R_temp)
            A_int=interp1(R_temp,A_temp,para.range_receiver, 'linear',nan);
            cnt=cnt+1;
        else % double eigenray
            clear A_int;
            peakind=find(diff(sign(diff(R_temp))))+1; peakind=[1,peakind.',length(R_temp)];
            for ind=1:length(peakind)-1
                A_int(ind,:)=interp1(R_temp(peakind(ind):peakind(ind+1)),A_temp(peakind(ind):peakind(ind+1)),para.range_receiver, 'linear',nan);
            end;
            cnt=cnt+length(peakind)-1;
        end;
    elseif length(A_temp)==1 && para.searchrays
        para.dense=0;
        if para.dense==1
            ang_foc=A_temp;
            iang=find(abs(para.scan-ang_foc)<5*eps);
            
            if iang > 1 && iang < length(para.scan)
                angs=(para.scan(iang-1:2:iang+1)-ang_foc)/2+ang_foc;
            elseif iang==1
                angs=(para.scan(iang+1)-ang_foc)/2*(-1:2:1)+ang_foc;
            elseif iang==length(para.scan)
                angs=(ang_foc-para.scan(iang-1))/2*(-1:2:1)+ang_foc;
            end;
            ok=0;loopcnt=0;
            while ~ok && loopcnt<=5
                newpara.scan=angs;
                [Eig, Count, sufbot, newRays] = tracerays(env, newpara);
                [T_temp R_temp A_temp TA_temp label] =selectarrival( categ_acs(inewklass,1),categ_acs(inewklass,2),categ_acs(inewklass,3),categ_acs(inewklass,4),categ_acs(inewklass,5), Eig,Count);
                if ~isempty(A_temp)
                    ok=1;
                    sel=find(any(label.'==1)); newK1=length(sel);  newK2=size(class,2)-size(label,2);
                    HIST.TIME      =[HIST.TIME;    Eig.TIME(sel,:),    nan*ones(newK1,newK2)];
                    HIST.RANGE     =[HIST.RANGE;   Eig.RANGE(sel,:),   nan*ones(newK1,newK2)];
                    HIST.THETA     =[HIST.THETA;   Eig.THETA(sel,:),   nan*ones(newK1,newK2)];
                    HIST.ANGLE     =[HIST.ANGLE;   Eig.ANGLE(sel,:),   nan*ones(newK1,newK2)];
                    HIST.DISTANCE  =[HIST.DISTANCE;Eig.DISTANCE(sel,:),nan*ones(newK1,newK2)];
                    COUNT.BOTTOM  =[COUNT.BOTTOM;Count.BOTTOM(sel,:),nan*ones(newK1,newK2)];
                    COUNT.SURFACE =[COUNT.SURFACE;Count.SURFACE(sel,:),nan*ones(newK1,newK2)];
                    COUNT.TURNING_ABOVE=[COUNT.TURNING_ABOVE;Count.TURNING_ABOVE(sel,:),nan*ones(newK1,newK2)];
                    COUNT.TURNING_BELOW=[COUNT.TURNING_BELOW;Count.TURNING_BELOW(sel,:),nan*ones(newK1,newK2)];
                    SUFBOT.BOTTOM_ANGLE =[SUFBOT.BOTTOM_ANGLE;sufbot.BOTTOM_ANGLE(sel,:),         1000*ones(newK1,size(SUFBOT.BOTTOM_ANGLE,2) -size(sufbot.BOTTOM_ANGLE,2))];
                    SUFBOT.BOTTOM_RANGE =[SUFBOT.BOTTOM_RANGE;sufbot.BOTTOM_RANGE(sel,:),  2*para.Rmax*ones(newK1,size(SUFBOT.BOTTOM_RANGE,2) -size(sufbot.BOTTOM_RANGE,2))];
                    SUFBOT.SURFACE_RANGE=[SUFBOT.SURFACE_RANGE;sufbot.SURFACE_RANGE(sel,:),2*para.Rmax*ones(newK1,size(SUFBOT.SURFACE_RANGE,2)-size(sufbot.SURFACE_RANGE,2))];
                    SUFBOT.SURFACE_ANGLE=[SUFBOT.SURFACE_ANGLE;sufbot.SURFACE_ANGLE(sel,:),       1000*ones(newK1,size(SUFBOT.SURFACE_ANGLE,2)-size(sufbot.SURFACE_ANGLE,2))];
                    
                    para_bak.scan=[para_bak.scan,angs(sel)];
                    [para_bak.scan,order]=sort(para_bak.scan);
                    HIST.TIME      = HIST.TIME(order,:);
                    HIST.RANGE     = HIST.RANGE(order,:);
                    HIST.THETA     = HIST.THETA(order,:);
                    HIST.ANGLE     = HIST.ANGLE(order,:);
                    HIST.DISTANCE  = HIST.DISTANCE(order,:);
                    COUNT.BOTTOM  = COUNT.BOTTOM(order,:);
                    COUNT.SURFACE = COUNT.SURFACE(order,:);
                    COUNT.TURNING_ABOVE=COUNT.TURNING_ABOVE(order,:);
                    COUNT.TURNING_BELOW=COUNT.TURNING_BELOW(order,:);
                    SUFBOT.BOTTOM_ANGLE = SUFBOT.BOTTOM_ANGLE(order,:);
                    SUFBOT.BOTTOM_RANGE = SUFBOT.BOTTOM_RANGE(order,:);
                    SUFBOT.SURFACE_RANGE= SUFBOT.SURFACE_RANGE(order,:);
                    SUFBOT.SURFACE_ANGLE= SUFBOT.SURFACE_ANGLE(order,:);
                    if ~isempty(Rays.X)
                        raylength=size(Rays.X,2);newraylength=size(newRays.X,2);
                        NEWX=nan(newK1,raylength);
                        NEWZ=nan(newK1,raylength);
                        NEWX(:,1:min(raylength,newraylength))=newRays.X(sel,1:min(raylength,newraylength));
                        NEWZ(:,1:min(raylength,newraylength))=newRays.Z(sel,1:min(raylength,newraylength));
                        Rays.X=[ Rays.X; NEWX ];
                        Rays.Z=[ Rays.Z; NEWZ ];
                        Rays.X            = Rays.X(order,:);
                        Rays.Z            = Rays.Z(order,:);
                    end;
                    
                    [~, R_temp, A_temp, ~, label] =selectarrival( categ_acs(inewklass,1),categ_acs(inewklass,2),categ_acs(inewklass,3),categ_acs(inewklass,4),categ_acs(inewklass,5), HIST,COUNT);
                    if issorted(R_temp) || issorted(-R_temp)
                        A_int=interp1(R_temp,A_temp,para.range_receiver, 'linear',nan);
                        cnt=cnt+1
                        
                    else % double eigenray
                        disp(['double'])
                        clear A_int;
                        peakind=find(diff(sign(diff(R_temp))))+1; peakind=[1,peakind.',length(R_temp)];
                        for ind=1:length(peakind)-1
                            A_int(ind,:)=interp1(R_temp(peakind(ind):peakind(ind+1)),A_temp(peakind(ind):peakind(ind+1)),para.range_receiver, 'linear',nan);
                        end;
                        cnt=cnt+length(peakind)-1;
                    end;
                else
                    angs=(angs-ang_foc)/2+ang_foc;
                    loopcnt=loopcnt+1;
                end;
            end;
            if ok
                class=[class;ones(newK1,1)*class(iang,:)];   class2=[class2;ones(newK1,1)*class2(iang,:)];
                class=class(order,:);   class2=class2(order,:);
            end;
        end;
    elseif ~para.searchrays
        A_int=[];
    end;
    if cntold<cnt
        startangle(:,cntold+1:cnt) = A_int.';
        [~,ir]=find(label);
        count(:,cntold+1:cnt)      = [categ_acs(inewklass,:),17+inewklass,1,ir(1)].'*ones(1,cnt-cntold);
    end;
end;

HIST.class  = class;
HIST.class2 = class2;
startangle = startangle(:,1:cnt);
HIST.count  = count(:,1:cnt);


%%
function[time,range,angle,angle_target,label]=selectarrival(bottom,surface,turning_above,turning_below,direction,HIST,COUNT)

eig_range=HIST.RANGE;
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
angle=eig_theta(label);angle=angle(~isnan(angle));
angle_target=eig_angle(label);angle_target=angle_target(~isnan(angle_target));
if direction>0;
    time=time(angle>0); range=range(angle>0); angle_target=angle_target(angle>0); angle=angle(angle>0);
    label(label>0 & eig_theta<0)=0;
end;
if direction<0;
    time=time(angle<0); range=range(angle<0); angle_target=angle_target(angle<0); angle=angle(angle<0);
    label(label>0 & eig_theta>0)=0;
end;

function segind = rayseg(X,range)
segind=zeros(1,length(range)+2);
segind(1)=1;
for iseg=1:length(range)
    [~,segind(iseg+1)]=min(abs(X-range(iseg)));
end;
segind(end)=length(X);

