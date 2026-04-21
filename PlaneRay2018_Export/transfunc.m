
%%   Authors: Jens M. Hovem and Shefeng Yan
%%   Revised: 2014/06/15
%%   Calculates the frequency domain function and store on disk

function [transfer_function,trans,c_red,eigen] = transfunc(HIST, COUNT, SUFBOT,env, para,origin)
% Calculates the transfer function by coherent additions of all eigenray contributions.
% The output varaiables are:
% transfer_function the total frequency domain transferfunction;
% c_red the applied reduction sound speed
% eigen= calculated by the function planerayfield and containing the applied values of
% startangle,timedelay,rayamp,distance,count, and range.
% The results are store on disk under in the struct "trans", containing
% trans.total= transfer_function;
% trans.direct = transfer_function_direct;
% trans.turn = transfer_function_turn;
% trans.surface=transfer_function_surface;
% trans.bottom=transfer_function_bottom;
% trans.sufbot=transfer_function_sufbot;
% trans.other=transfer_function_other;
% origin==1 for the calculation of the transmission loss
% origin==2 for the calculation of the pulseresponse


BOTTOM_ANGLE=SUFBOT.BOTTOM_ANGLE;
BOTTOM_RANGE=SUFBOT.BOTTOM_RANGE;
SURFACE_ANGLE=SUFBOT.SURFACE_ANGLE;
SURFACE_RANGE=SUFBOT.SURFACE_RANGE;

if ~para.quiet,disp('Start calculating transfer function ...');tic;end;
eigen = planerayfield(HIST, COUNT, env,para);
startangle=eigen.startangle;
timedelay=eigen.timedelay;
rayamp=eigen.rayamp;
distance=eigen.distance;
count=eigen.count;
phone_positions=para.range_receiver;

N_phones=length(phone_positions);
R_bathy = env.bathy(1,:);
Z_bathy = env.bathy(2,:);
z=env.ssp(:,1);
c=env.ssp(:,2);
del_z=z(2)-z(1);
start_angle=para.start_angles;
bottom_depth=Z_bathy(1);
c_depth=c(fix(bottom_depth/del_z)-1);

c_red=para.c_red;


sigma_surface = env.sigma_surface; %Sea surface roughness
sigma_bottom = env.sigma_bottom; %Bottom roughness
konst=8.686*2*pi;
dB_lambda_p1=env.ap1;  alpha_p1=dB_lambda_p1/konst;
dB_lambda_p2=env.ap2;  alpha_p2=dB_lambda_p2/konst;
dB_lambda_s2=env.as2;  alpha_s2=dB_lambda_s2/konst;
ccp1=env.cp1*(1+1i*alpha_p1);
ccp2=env.cp2*(1+1i*alpha_p2);
ccs2=env.cs2*(1+1i*alpha_s2);
rrho1=env.rho1;
rrho2=env.rho2;
RR=env.RR;
sedthick=env.laythick;
% end geoacoustic model

f=para.frequency;
cf=para.carrier_frequency;
frequency=f+cf;
nfreq=length(f);
omega=2*pi*f;
 if length(RR)<=2 % range independent
     cp1=ccp1(1); cp2=ccp2(1); cs2=ccs2(1);  rho1=rrho1(1); rho2=rrho2(1); d=sedthick(1);
 end;
transfer_function = zeros(N_phones,nfreq);
transfer_function_inc = zeros(N_phones,nfreq);
transfer_function_direct = zeros(N_phones,nfreq);
transfer_function_turn = zeros(N_phones,nfreq);
transfer_function_surface = zeros(N_phones,nfreq);
transfer_function_bottom = zeros(N_phones,nfreq);
transfer_function_sufbot = zeros(N_phones,nfreq);
transfer_function_other = zeros(N_phones,nfreq);
beam_spectrum=ones(1, nfreq);
T=env.temperature;
S=env.salinity;
ph=8;% fixed ph value)
D=mean([ para.z_receiver para.z_source]);

[temp, a a1 a2 a3] = frangar(T,S,D,ph,(cf+f)/1000);
alpha=a/1000;

lambda=1500./frequency;
for jj=1:N_phones

    
    receiver_range=phone_positions(jj);
    for cnt=1:size(startangle,2);
        current_angle=startangle(jj,cnt);
        theta=current_angle*pi/180;
        if origin==2;
            beam_spectrum  = getbeamspec(para, theta,1);
        end;
 
        if ~isnan(current_angle)
            categ=count(:,cnt);
            N_bottom  = categ(1); N_surface = categ(2);
            N_above   = categ(3); N_below   = categ(4);
            if N_bottom>para.bothit,break;end;
            surface_ref=1; bottom_ref=1; aboveturn_ref=1; belowturn_ref=1;
            amp=rayamp(jj,cnt);
            
            attenuation=10.^(-distance(jj,cnt)*alpha/20);
            delay=timedelay(jj,cnt)-(receiver_range- para. range_source)/c_red;
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
                else
                    for pp=1:N_bottom
                        bottom_range=range_bottom(pp);
                        theta=theta_bottom(pp);
                        if length(R_bathy)>2 % not flat bottom
                            bottom_depth=Z_bathy(find(bottom_range-R_bathy>=0,1,'last'));
                            nn=(bottom_depth/del_z);
                            if  ~isempty(nn);   c_depth=c(fix(nn));end
                        end;
                        if length(RR)>2 % range dependent
                            ind = find(bottom_range-RR>=0, 1,'last');
                            cp1=ccp1(ind); cp2=ccp2(ind); cs2=ccs2(ind);  rho1=rrho1(ind); rho2=rrho2(ind); d=sedthick(ind);
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
            trf = E.*amp.*surface_ref.*bottom_ref.*aboveturn_ref.*belowturn_ref.*attenuation.*beam_spectrum;
          trf_inc=abs(trf);
            jjj=isnan(theta);
            
            if jjj>0; disp([ 'help  '  num2str(jj) '    ' num2str(receiver_range) '    '  num2str(theta) ]);end;
            transfer_function(jj,:)=transfer_function(jj,:) +  trf;
             transfer_function_inc(jj,:)=transfer_function_inc(jj,:) +  trf_inc;
            switch categ(6)
                case {1}
                    transfer_function_direct(jj,:)=transfer_function_direct(jj,:) +  trf;
                case {2,3,4,5}
                    transfer_function_turn(jj,:)=transfer_function_turn(jj,:) +  trf;
                case {6,7,8,9}
                    transfer_function_sufbot(jj,:)=transfer_function_sufbot(jj,:) +  trf;
                case {10,11,12,13}
                    transfer_function_bottom(jj,:)=transfer_function_bottom(jj,:) +  trf;
                case {14,15,16,17}
                    transfer_function_surface(jj,:)=transfer_function_surface(jj,:) +  trf;
                    transfer_function_other(jj,:)=transfer_function_other(jj,:) +  trf;
            end;
        end;
    end;
end;% jj
trans.total= transfer_function;
trans.direct = transfer_function_direct;
trans.turn = transfer_function_turn;
trans.surface=transfer_function_surface;
trans.bottom=transfer_function_bottom;
trans.sufbot=transfer_function_sufbot;
trans.other=transfer_function_other;
trans.inc=transfer_function_inc;
save trans trans
if ~para.quiet,disp(['End calculating transfer function. Total time: ',num2str(toc),' s']);end;


%


