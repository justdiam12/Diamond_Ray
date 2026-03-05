clf
clear all
load ssp
% load SourcePositions
RR = 5460;
SD = 55;
zw = ssp_depths;
cw = ssp;
% cw = Cssp;
dcz = spline((zw(1:end-1)+zw(2:end))/2,diff(cw),zw);
cb = 1610;
zb = 0;
subplot(151)
plot(cw,zw);axis ij
% ylim([0 20])
xlabel('sound speed (m/s)')
%Manual SSP
% zw = linspace(0,80,100);
% cw = zw*0+1500;
% dcz = spline((zw(1:end-1)+zw(2:end))/2,diff(cw),zw);
% cb = 1610;
% zb = 0;
% 
subplot(1,5,[2:5])
RD = 55;
H = 200;
zw(end)=H;

%DP initial guess (ISO)
lnchang = atand((SD-RD)./RR);
thtDP = lnchang;
%SP initial guess (ISO)
lnchang = atand((SD+RD)./RR);
thtSP = lnchang;
%BP initial guess (ISO)
lnchang = atand((SD-(2*H-RD))./RR);
thtBP = lnchang;
%BSP initial guess (ISO)
lnchang = atand((SD-(2*H+RD))./RR);
thtBSP = lnchang;
%SBP initial guess (ISO)
lnchang = atand((SD+(2*H-RD))./RR);
thtSBP = lnchang;
%BSBP initial guess (ISO)
lnchang = atand((SD-(4*H-RD))./RR);
thtBSBP = lnchang;
%SBSP initial guess (ISO)
lnchang = atand((SD+(2*H+RD))./RR);
thtSBSP = lnchang;
% 
[Range,Depth,time] = eigenray_waveguide([thtDP thtSP thtBP thtBSP thtSBP-.5]-2,...
                     .2,.1,RR,RD,SD,dcz,cw,zw,cb,zb);
clf
plot([0 RR],[H H],'k',[0 RR],[0 0],'b');axis ij
hold on             
x(1) = plot(RR-Range(:,1),Depth(:,1),'k-','Linewidth',6);
x(2) = plot(RR-Range(:,2),Depth(:,2),'k-','Linewidth',4);
x(3) = plot(RR-Range(:,3),Depth(:,3),'k:','Linewidth',4);
x(4) = plot(RR-Range(:,4),Depth(:,4),'k-','Linewidth',2);
x(5) = plot(RR-Range(:,5),Depth(:,5),'k:','Linewidth',2);
% theta =-[0.21862 -1.2565 3.9271 4.8238 -4.8225];
% [Range,Depth,time] = rayfan_waveguide(theta,1,RR,SD,dcz,cw,zw,cb,zb);
% plot(Range,Depth);axis ij
