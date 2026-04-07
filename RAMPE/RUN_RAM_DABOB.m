%%%LOAD ENVIRONMENTAL PARAMETERS%%%%
load("DABOBathy.mat");
[dl] = -interp2(lon,lat,d,linspace(-122-48.175/60,-122.842 ,5000),linspace(47+46.378/60, 47.73,5000));
rl = distance(47+46.378/60,-122-48.175/60,linspace(47+46.378/60, 47.73,5000),linspace(-122-48.175/60,-122.842 ,5000),wgs84Ellipsoid);
r = [0:1:rl(end)];
d = interp1(rl,dl,r);

load("ssp.mat");
cw = ssp;
zw = ssp_depths;

freqs = [3450:3550];
for freqdex = [1:length(freqs)];
%%%%% place into format for RAM-PE %%%%
fprintf(int2str(freqs(freqdex)))
RNG =  [r];
DPTH = [d]; 

%%
c0 = 1500;  %Reference Sound Speed
freq = freqs(freqdex); %Source Frequency
SD = DPTH(1)-4;     %Source Depth
RD = 5;     %Receiver Depth for PEoutput.line
RR = RNG(end);    %Receiver Range
dz =  .05;    %Depth Grid Spacing
dr =  .05;    %Range Grid Spacing 
ndr = 20;     %Number of Range Outputs
ndz = 20;     %Number of Depth Outputs

cw = cw;
zw = zw;

attnw = cw*0;
rhow = cw*0+1;

zb = zw(end)+[0 400 500 600 700];%Depth Bottom
cb =1700*ones(size(zb)); %Sound Speed Bottom

rhob = [1 1 1 1 1]*1.50;           %Density Bottom
attnb = 0.2 * [1 1 2 5 100]; %in dB/lambda
zmplt = 600;                 %Maximum Depth of Grid Output


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adding AIR layer - Use the form Row Vector ca = [... 340 340 340 340]
%                                           za = [... za3 za2 za1 0]
za =   -[150 20 15 10 dz];
ca =    [340 340 340 340 340];
rhoa =  [0.00128 0.00128 0.00128 0.00128 0.00128];
attna = [100 10 1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load in Bathymetry
%%%%%%%%%%%%%%%%%%%%%%%%%%
rb=RNG;
z_rb=DPTH-za(1)+dz;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load in Rough Surface
%%%%%%%%%%%%%%%%%%%%%%%%%%
x = RNG(1):100:RNG(end);
surface = x*0+sin(x*2*pi/1000)*0;

CppInputWriterRS(freq,SD,RD,RR,dr,dz,ndr,ndz,rb,z_rb,zw,cw,rhow,attnw,za,ca,rhoa,attna,zb,cb,rhob,attnb,zmplt,c0,surface,x);


%%%%%%%%%%%%%%%
% RUN PE CODE %
%%%%%%%%%%%%%%%
fprintf('Running FORTRAN CODE...')
tic
status = dos('ramrun.exe');
toc
% readPEsstDRD
readPEsstVECTOR
imagesc(Range,Depth+za(1),20*log10(abs(P)))
ylim([-20 200]);caxis([-80 -40])
drawnow

z = Depth+za(1);
zdex = find(Depth>-10&Depth<200);
Psv(:,:,freqdex)=P(zdex,:);

end
