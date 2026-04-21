function plotenv(env,para,movedown)
%   Authors: Jens M. Hovem and Shefeng Yan
%   Revision: July 2013

if nargin<3
    movedown=0;
end;

ssp          = env.ssp;
R_bathy      = env.bathy(1,:);
Z_bathy      = env.bathy(2,:);
R_slope      = env.bathslope(1,:);
laythick     = env.laythick;
R_max        = para.Rmax;
z_source     = para.z_source;
z_receiver   = para.z_receiver;
range_source = para.range_source;
RR           = env.RR;
RR           = [RR(1:length(laythick)),R_max];

z=ssp(1:end-1,1);
c=ssp(1:end-1,2);
c_min=5*fix(min(c)/5)-5;
c_max=5*fix(max(c)/5)+10;
z_bottom=max(z);
theta_min=min(para.start_angles);
theta_max=max(para.start_angles);
if ~isfield(para,'plotssp')
    para.plotssp=1;
end;
if para.plotssp==1;
    
    fonts=14; %lines=2;
    set(gcf,'DefaultTextFontSize', fonts,'DefaultAxesFontSize', fonts);
    figposition =get(gcf,'position');
    figposition(1)=figposition(1)-figposition(3)/4;
    figposition(3)=figposition(3)+figposition(3)/2;
    figposition(2)=figposition(2)-figposition(4)*movedown;
    figposition(4)=figposition(4)+figposition(4)*movedown;
    
    
    rect1=[0.1 0.2 0.15 0.4];
    rect2=[0.3 0.2 0.6 0.4];
    set(gcf,'position',figposition);
    axes('position',rect1); % Sound speed profile
    plot(c,z,'linewidth',1);
    xlabel('Sound speed – m/s')
    ylabel('Depth – m');
    axis([c_min, c_max,  0, z_bottom]);
    set(gca,'ydir','reverse');
    title(env.title);
    
    axes('position',rect2); % plotting of bottom topography and layering
end;
laythick=[laythick,laythick(end)];
d=zeros(1,length(R_slope));
for ir=1:length(R_slope)
    btrange=R_slope(ir);
    ind = find(btrange-RR>=0, 1,'last');
    d(ir)=laythick(ind);
end;

newRR=[R_slope(2:end-1);R_slope(2:end-1)];newRR=newRR(:).';
newRR=[R_slope(1),newRR,R_slope(end)];
newthick=[d(2:end-1);d(2:end-1)];newthick=newthick(:).';
newthick=[d(1),newthick,d(end)];
newZZ=depth(newRR,R_bathy,Z_bathy)+newthick;
z_bottom=max(max(newZZ),z_bottom);

X_fill_1=[0,R_bathy,R_max,0];
Y_fill_1=[z_bottom,Z_bathy,z_bottom,z_bottom];
X_fill_2=[0,newRR,R_max 0];
Y_fill_2=[z_bottom,newZZ, z_bottom, z_bottom];
frame_x=[0,R_max,R_max,0,0];
frame_y=[z_bottom,z_bottom,0,0,z_bottom];
rockcolor=[1.0 0.69 0.39];
sedimentcolor=[0.95 0.87 0.73];
fill(X_fill_1/1000, Y_fill_1, sedimentcolor);hold on;

if any(laythick~=0)
    fill(X_fill_2/1000, Y_fill_2, rockcolor);
    
else
    fill(X_fill_1/1000, Y_fill_1, rockcolor);
end;
lines=2;
plot(R_bathy/1000, Z_bathy,'k',frame_x/1000,frame_y,'k',range_source/1000, z_source,'rp', 'linewidth',lines);
xlabel('Range – km');
axis([0,R_max/1000,0,z_bottom]);
set(gca,'ydir','reverse');
if para.plotssp==1
    set(gca,'YAxisLocation','right');
end;
title(['Source depth = ' num2str(z_source) ' m, ',  '  Angles= ',num2str(fix(theta_min)) '{\circ}',' : ',num2str(fix(theta_max)),'{\circ}']);


function depth_current = depth(range_current, R_bathy, Z_bathy)
depth_current=zeros(1,length(range_current));
for ii=1:length(range_current)
    depth_current(ii)=Z_bathy(find(range_current(ii)-R_bathy>=0,1,'last'));
end;