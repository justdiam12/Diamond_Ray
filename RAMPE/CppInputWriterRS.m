%Write input file for a RamPE_RS (WVGDinput.in)

%Notes: Compile/Build RAM_PE.an with correct input filename

function CppInputWriterRS(freq,SD,RD,RR,dr,dz,ndr,ndz,rb,z_rb,zw,cw,rhow,attnw,za,ca,rhoa,attna,zb,cb,rhob,attnb,zmplt,c0,surface,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writing Waveguide Properties input file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['ram.in'];
% filename = [pwd '/WVGDinput.in'];
fid = fopen(filename,'wt');

fs = freq;
zs = SD-za(1);
zr = zs;
rmax = RR;
zmax = zb(end)-za(1);
np = 8;
ns = 1;
rs = 0;%1*max([ca cw cb])./freq;
k0 = 2*pi*fs/c0;

isrc_type = 0;  %0 IS MONOPOLE, 1 IS BEAM
r0  = 0.0; 
r1  = 0.0;  
theta0 = 5;              
beam_width = 15%100*1/k0;
%Writing Header
fprintf(fid,['Test input file for range-dependent RAM_PE' '\n']);
fprintf(fid,'%d %6.6f %6.6f %6.6f %6.6f\n',[isrc_type r0 r1 theta0 beam_width]);
fprintf(fid,'%6.6f %6.6f %6.6f\n',[fs zs zr]);
fprintf(fid,'%6.6f %6.6f %6.0f\n',[rmax dr ndr]);
fprintf(fid,'%6.6f %6.6f %6.0f %6.6f\n',[zmax dz ndz zmplt]);
fprintf(fid,'%6.6f %6.0f %6.0f %6.0f\n',[c0 np ns rs]);

%Bathemetry Information
%% INPUT
%Range Vector:      rb
%Bottom Depth:      zb

RINDP = [rb;z_rb];
fprintf(fid,'%6.6f %6.6f\n',RINDP);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%Surface Loop
%save water column properties
zasave = za;
casave = ca;
rhoasave = rhoa;
attnasave = attna;
zwsave = zw;
cwsave = cw;
rhowsave = rhow;
attnwsave = attnw;

for hdex = 1:length(x) %Number of surface height definitions in PE computational Range

za = zasave;
ca = casave;
rhoa = rhoasave;
attna = attnasave;
zw = zwsave; %Reload Water Column Properties
cw = cwsave;
rhow = rhowsave;    
attnw = attnwsave;

h = -surface(hdex);
za(end) = za(end)-h;
za(end+1) = za(end)+dz; %Not sure if this dz is needed, check discontinuity (helps smoove)
ca(end+1) = cw(1);
rhoa(end+1) = rhow(1);
attna(end+1) = attnw(1);

zdex = find(zw<=za(end));

    zw(zdex)=[];
    cw(zdex)=[];
    rhow(zdex)=[];
    attnw(zdex)=[];

    %AIR/Water Column Profile
%% INPUT (Currently set for Pekeris Waveguide Case)
%Depth Vector:      za,zw
%Sound Speed:       ca,cw
%Density:           rhoa,rhow
%Attn:              attna,attnw
if hdex>1
    fprintf(fid,'%7.2f             rp\n', x(hdex));
end
%AIR/Water Soundspeed Profile
SSPFL_AW = [[za zw]-za(1);[ca cw]];
fprintf(fid,'%6.6f %6.6f\n',SSPFL_AW);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%AIR/Water Density Profile
DENPFL_AW = [[za zw]-za(1);[rhoa rhow]];
fprintf(fid,'%6.6f %6.6f\n',DENPFL_AW);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%AIR/Water Attenuation Profile
ATTNPFL_AW = [[za zw]-za(1);[attna attnw]];
fprintf(fid,'%6.6f %6.6f\n',ATTNPFL_AW);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%Seabed Profile
%% INPUT (Currently set for Pekeris Waveguide Case)
%Depth Vector:      zb
%Sound Speed:       cb
%Density:           rhob
%Attn:              attnb

%Seabed Soundspeed Profile
SSPFL_B = [[zb]-za(1);cb];
fprintf(fid,'%6.6f %6.6f\n',SSPFL_B);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%Seabed Density Profile
DENPFL_B = [[zb]-za(1);rhob];
fprintf(fid,'%6.6f %6.6f\n',DENPFL_B);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

%Seabed Attenuation Profile
ATTNPFL_B = [[zb]-za(1);attnb];
fprintf(fid,'%6.6f %6.6f\n',ATTNPFL_B);
fprintf(fid,'%6.0f %6.0f\n',[-1 -1]);

if hdex==1; %PRINT GHOST in FILE FOR REFERENCE
filename2 = ['WVGDinput_ghost.in'];
fid2 = fopen(filename2,'wt');
fprintf(fid2,['Test input file for range-dependent RAM_PE' '\n']);
fprintf(fid2,'%6.6f %6.6f %6.6f\n',[fs zs zr]);
fprintf(fid2,'%6.6f %6.6f %6.0f\n',[rmax dr ndr]);
fprintf(fid2,'%6.6f %6.6f %6.0f %6.6f\n',[zmax dz ndz zmplt]);
fprintf(fid2,'%6.6f %6.0f %6.0f %6.0f\n',[c0 np ns rs]);
fprintf(fid2,'%6.6f %6.6f\n',RINDP);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',SSPFL_AW);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',DENPFL_AW);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',ATTNPFL_AW);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',SSPFL_B);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',DENPFL_B);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fprintf(fid2,'%6.6f %6.6f\n',ATTNPFL_B);
fprintf(fid2,'%6.0f %6.0f\n',[-1 -1]);
fclose(fid2);
end
end

fclose('all');

