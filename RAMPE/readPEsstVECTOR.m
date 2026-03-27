%readPEsst

%Parse the input file for parameters of output grid
input_filename='ram.in';
%  be sure above matches the RAM input file
fid_input=fopen(input_filename,'r');
tline = fgetl(fid_input); %reads the line as entire string
tline = fgetl(fid_input); %reads the line as entire string
tline = fgetl(fid_input); %first string containing input data
%[val]=strread(tline,'%f',2);freq=val(1);zs=val(2); %pars out the string
%change 14MAY
[val]=strread(tline,'%f',3);freq=val(1);zs=val(2);zr=val(3); %pars out the string
tline = fgetl(fid_input);
[val]=strread(tline,'%f',3);rmax=val(1);dr=val(2);ndr=val(3);
tline = fgetl(fid_input);
[val]=strread(tline,'%f',4);zmax=val(1);dz=val(2);ndz=val(3);zmplt=val(4);
tline = fgetl(fid_input);
[val]=strread(tline,'%f',4);c0=val(1);np=val(2);ns=val(3);rs=val(4);
tline = fgetl(fid_input);
[val]=strread(tline,'%f',2);rbs=val(1);zbs=val(2);
fclose(fid_input);
% convert to operational variables (i.e., what we record & plot)
dr=dr; 
rge=rmax;
dz=dz;
dep=zmplt;
k0=2*pi*freq/c0;

%  what the input file should look like....
%range-dependent example
%300.0 1.0 42.0         freq zs zr
%5000.0 5.0 1          rmax dr ndr
%500.0 1 1 150.0       zmax dz ndz zmplt
%1532.5 8 1 0.0           c0 np ns rs
%0.0 75.0                rb zb


%%%%%%%%%%%MOD - David Dall'Osto 12/17/2018%%%%%%%%%%%%%%%%%
%BUILD RANGE VECTOR AND DEPTH VECTOR UTILIZED BY THE PE build in SST
%    UPDATED FOR 4 FILE INPUT 4/23/2019                    %
rhow = 1000; %assumption of constant impedence here

dat_filename='tl.grid'; 
fid = fopen('tl.grid','r');
PEdata1 = fread(fid,inf,'float32');
fclose(fid)
PEdata = PEdata1;
wrapdex = find(PEdata == PEdata(4));
PEdata(wrapdex) = [];
PEdata(1:3)=[];
xsize = (diff(wrapdex(1:2))-1);
ysize = length(PEdata)/xsize/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i
OUT=reshape(PEdata,xsize,2,ysize);
OUT(:,:,1)=[]; %FIRST IS REPEAT
OUT = squeeze(OUT(:,1,:)+i*OUT(:,2,:));
Depth=[dz:ndz*dz:xsize*ndz*dz]-dz; %lost points at surface and bottom
Range=[ndr*dr:ndr*dr:dr*ndr*(ysize-1)];

%reallocate range e(ik0r) phase dependence
P = OUT.*repmat(exp(i*k0*Range-i*pi/4),length(Depth),1);


dat_filename='tl.up'; 
fid = fopen('tl.up','r');
PEdata2 = fread(fid,inf,'float32');
fclose(fid)
PEdata = PEdata2;
PEdata(wrapdex-3) = [];
xsize = (diff(wrapdex(1:2))-1);
ysize = length(PEdata)/xsize/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i
OUT=reshape(PEdata,xsize,2,ysize);
OUT(:,:,1)=[]; %FIRST IS REPEAT
OUT = squeeze(OUT(:,1,:)+i*OUT(:,2,:));
Depth=[dz:ndz*dz:xsize*ndz*dz]-dz; %lost points at surface and bottom
Range=[ndr*dr:ndr*dr:dr*ndr*(ysize-1)];

%reallocate range e(ik0r) phase dependence
Pm = OUT.*repmat(exp(i*k0*Range-i*pi/4),length(Depth),1);


dat_filename='tl.down'; 
fid = fopen('tl.down','r');
PEdata3 = fread(fid,inf,'float32');
fclose(fid)
PEdata = PEdata3;
%wrapdex = find(PEdata == PEdata(4));
PEdata(wrapdex-3) = [];
xsize = (diff(wrapdex(1:2))-1);
ysize = length(PEdata)/xsize/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i
OUT=reshape(PEdata,xsize,2,ysize);
OUT(:,:,1)=[]; %FIRST IS REPEAT
OUT = squeeze(OUT(:,1,:)+i*OUT(:,2,:));
Depth=[dz:ndz*dz:xsize*ndz*dz]-dz; %lost points at surface and bottom
Range=[ndr*dr:ndr*dr:dr*ndr*(ysize-1)];

%reallocate range e(ik0r) phase dependence
Pp = OUT.*repmat(exp(i*k0*Range-i*pi/4),length(Depth),1);

% 

dat_filename='tl.left'; 
fid = fopen('tl.left','r');
PEdata4 = fread(fid,inf,'float32');
fclose(fid)
PEdata = PEdata4;
wrapdex = find(PEdata == PEdata(1));
PEdata(wrapdex) = [];
xsize = (diff(wrapdex(1:2))-1);
ysize = length(PEdata)/xsize/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i
OUT=reshape(PEdata,xsize,2,ysize);
OUT(:,:,1)=[]; %FIRST IS REPEAT
OUT = squeeze(OUT(:,1,:)+i*OUT(:,2,:));
Depth=[dz:ndz*dz:xsize*ndz*dz]-dz; %lost points at surface and bottom
RangeL=[ndr*dr:ndr*dr:dr*ndr*(ysize-1)]-dr;

%reallocate range e(ik0r) phase dependence
Pl = OUT.*repmat(exp(i*k0*RangeL-i*pi/4),length(Depth),1);


if prod(size(Pl))~=prod(size(P))
Pl(:,end)=[];
RangeL(end)=[];
end


% 
% 
dat_filename='tl.right'; 
fid = fopen('tl.right','r');
PEdata5 = fread(fid,inf,'float32');
fclose(fid)
PEdata = PEdata5;
wrapdex = find(PEdata == PEdata(1));
% PEdata(1:(xsize+2)*2) = [];
PEdata(wrapdex) = [];
xsize = (diff(wrapdex(1:2))-1);
ysize = length(PEdata)/xsize/2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%i
OUT=reshape(PEdata,xsize,2,ysize);
OUT(:,:,1)=[]; %FIRST IS REPEAT
OUT = squeeze(OUT(:,1,:)+i*OUT(:,2,:));

OUT = circshift(OUT,[0 -1]); %what's this? seem's neccessary in build as of 5/3/2019

Depth=[dz:ndz*dz:xsize*ndz*dz]-dz; %lost points at surface and bottom
RangeR=[ndr*dr:ndr*dr:dr*ndr*(ysize-1)]+dr;

%reallocate range e(ik0r) phase dependence
Pr = OUT.*repmat(exp(i*k0*RangeR-i*pi/4),length(Depth),1);
% 

% 
if prod(size(Pr))~=prod(size(P))
Pr(:,end)=[];
RangeR(end)=[];
end
% 
% plot(Range,20*log10(abs(P(40,:))),'o',...
%      Range,20*log10(abs(Pp(40,:))),'.',...
%      Range,20*log10(abs(Pm(40,:))),'x',...
%      RangeL,20*log10(abs(Pl(40,:))),'s',...
%      RangeR,20*log10(abs(Pr(40,:))),'p','linestyle','none')
% 
% %  
 omega = k0*c0;
% %NOW CONSTRUCT VELOCITY COMPONENTS
dPr = (Pr-Pl)./repmat(RangeR-RangeL,length(Depth),1);
dPz = (Pp-Pm)/(2*dz);
Vr = -dPr./(-i*omega*rhow);
Vz = -dPz./(-i*omega*rhow);


Iz = real(P.*conj(Vz));
Ir = real(P.*conj(Vr));

Qz = imag(P.*conj(Vz));
Qr = imag(P.*conj(Vr));

% pcolor(Range,Depth,real(Ir));shading flat;caxis([-1 1]/1e12)
