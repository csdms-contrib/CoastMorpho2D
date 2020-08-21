clear; close all;clc
%NOTES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qs is the total lateral transport, E is the vertical erosion flux

%lateral boundary options
% 0 is closed (no flux if nothign specified. no-gradient if AW equal to 1 or -1)
% 1 is periodic

%options for lateral wave condtions
%AW

%A
%0: not in the domain
%1: a normal cell
%2: the open sea boundary
%10: the river boundary
%%%FALSE%3 and -3: no-gradient boundaries (for the lateral). 3 is left (1)  -3 is right (end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%to get always the same random numbers
rng(2); %can use any number in here


%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=struct;
P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.5;%0.1;%0.5;%0.2;%0.2;%0.1; % minimum water depth [m]. 0.2 is good NEEDS TO BE SMALLER THAN hwSea_lim
P.DiffSmud=1; %coefficient for tidal dispersion [-]. DO NOT CHANGE
P.DiffSsand=1;
P.DoMUD=1;%10;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)

%Sea level rise
P.RSLR=20/1000/365;  %from mm/yr to m/day (the time unit is the day!)
%P.RSLR=1.5/1000/365;  %from mm/yr to m/day (the time unit is the day!)


%Tide
P.Ttide=12.5/24; %tidal period [day]
P.Trange=1.5; %mean tidal Trange [m]
P.TrangeVEG=P.Trange;%tidal Trange for vegetation [m]. Generally same of tidal range

%Storm surge
P.alpha_surge=0.25;%1;%2*0.25;%0.25;%0.25;%how much the surge contributes to tidal prsim   0.1;%

%Swell waves
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
P.Pswelldir=0.5;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0; %if zero, only low angle waves
P.Ho=1.5;     %3;%1.9;%1.9;%1.7;%2; %boundary swell height Hs [m]
P.Tp_swell=8;%8;%8.2;%8;%6;% %boundary swell period Tp [m]
P.nrefrac=4;%either 0,1,2,3,4  Wave refraction. If zero there is no wave refraction
P.multifrequency=0;%on/off
P.wavediffraction=1;%on/off
P.Cbr=0.55;
P.Cbed=0.038;%wave bed friction if you use Jonswap (0.067 or 0.038)
P.wavefrictionCollins=0;

%Wind for sea waves
P.wind=7;%reference wind speed [m/s]

%Edge erosion
P.aw=0.3/365; %wave edge erodability m/yr/W/m2
P.maxedgeheight=2;
P.fox=0.5;%fraction of edge eroded material that is oxidized.

%Wind waves and swell numerics
P.hwSwell_lim=0.2;%0.10001;%2;%0.2 %limiter water depth for swell
P.hwSea_lim=0.2;%0.10001;%.2;%0.5;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE" ; NEEDS TO BE LARGER THAN KO!!!!

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=20/1000; %Sea boundary SSC for mud [g/l]
P.co3=0/1000; %Sea boundary SSC for mud [g/l]

%Impose the sediment discharge input at the river mouth
P.Qmouth=5; %river discharge per unit of cell [m2/s]
P.hmouth=5; %water depth [m]
%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
P.co2mouth=0;%500/1000; %SSC of mud at the river [g/l]

%Manning coeffinent unvegeated (same for sand and mud)
P.Cb=0.02;

%Downslope paramters for sand and mud (proportional to the sediment transport Qs!!!!)
P.alphaSAND=5;%2;%10;%5;%5;%3;%5;%3;%3;%10;%2; %coefficient for bedload downslope of sand. Calibrated with JMSE 2018, do not change
P.hlimC=1;%limit to apply to downslope of sand (and maybe also mud) Calibrated with JMSE 2018, do not change limiter height for total downlsope flux, to increase downslope in very shallow areas
P.downslopeSANDseawaves=10;% multiplication for sea-waves dowbsloep transport (^5 as for the swells) for SAND

%Sand
P.d50_1=0.25/1000/1;% %sand grain size [m]
P.ws1=0.02/1;% (note that sand with D50=500um has 0.05 %m/s)
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);

%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%settling velocity
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);

%Mud parameters
P.me=0.1*10^-4*24*3600;  %per day!!!
P.taucr=0.2;
P.tcrgradeint=0;% Pa/m
P.leveltauincrease=P.TrangeVEG/2;%1;
P.crMARSH=0.1/365;%creep coefficient vegetated
P.crMUD=3.65/365;%creep coeffcinet
P.alphaMUD=0.25; %coefficient for bedload downslope of mud. added April 2019. Similar to P.alphaSAND. Changed to 0.25 from 0.5 because initially used bulk1 instead of bulk2

%Correction for proceeses duration (to scale waves and tidal transport, because waves do not occur all the time!)
P.fMFswell=1; %use 1 if you use the equvilanet wave height (reduce it if you use the dominant wave height)
P.fMFsea=1; %use 1 if you use the equvilanet wind speed (reduce it if you use the dominant wind speed)
P.fMFriver=10/365;

%Vegetation parameters
P.dBlo=0;%lower limit veg
P.dBup=P.TrangeVEG/2;%upper limit veg
P.Cv=0.1;%Manning for vegetated ara
P.wsB=1/1000;%Mud Settling velocity for vegetated ara
P.taucrVEG=0.5;%Critical sheak stress for vegetated areas

%Organic accretion by vegetation
P.AccreteOrganic=1;
P.Korg=6/1000/365;%5/1000/365;%P.Korg=org/365;%5/1000/365;

%ON/OFF processes
P.VEGETATION=1;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.computemud=1;
P.computesand=1;
P.computeSwellwave=1;
P.computeSeaWaves=1;
P.computeEdgeErosionSwell=0;
P.computeEdgeErosionSea=1;
P.compute_currentbankerosion=0;
P.computetide=1;
P.computeriver=0;

%Correction for second-order river dynamics
P.riverwaterlevel=0;
P.rivermomemntumcorrection=0;

%Various boundary condtions
P.periodic=1;
P.imposeseaboundarydepthmorphoALL=0;%0; %to use when a channel mouth is at a boundary

%Ebb-flood momentum correction
P.ebbfloodcorrection=1;
P.residualcurrents=0;

%Pond dynamics
P.calculateponddynamics=0;
P.parameterforpond1=1;
P.parameterforpond2=1;
P.parameterforpond3=1;

%Stratigraphy
P.evolvestratigraphy=1;
P.VEGstratigraphy=0;%if 0 then you put the organic into the mud. If 1 then you calculate the organic as a sediment per se (advection, divergence,etc)
P.VEGonsand=0;

%Stratigraphy parameters
P.nlyr=20; %max number of layers
P.dlyr=1; %thickenss of layers
P.tlyrU=3; %max depth to add layer %must be larger than dlyr
P.tlyrD=1; %min depth merge layers %mus be smaller than dlyr
P.tcko=10;%10;%tickness of bed layer
P.levo=7;%intial level occupied
P.YUi=2;%initial thickess of active layer
P.initialfU=1;%initial composition of the active layer
P.initialf=1;%initial composition of all the layers

P.reducefractionsediment=1;%this should be 1 unless you to strange stuff with straigraphy

%Global numerical limiters
P.limitdeltaz=5;%m
P.limitmaxup=1;%m

%Time parameters
tmax=20000;%
tINT=1;%how many time steps you want to do the plot (if 1 you plot every time step). Does not affect the computation




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TIME SERIES

%1-Storm surges
numberserie=20000;%if you change this you will change the actual values in the time series, rememebr!
numberevents=numberserie/2;

lag=(2*365-20)*ones(numberevents,1);%exprnd(1.8*365,numberevents,1);%lag(lag<1)=1;
duration=20*ones(numberevents,1);%exprnd(2,numberevents,1); %duration=exprnd(0.26,numberevents,1);lag(lag<0.01)=0.01; %duration=exprnd(0.0003*365,numberevents,1);lag(lag<0.01)=0.01;
dtOserie=ones(numberserie,1);
dtOserie(1:2:end)=lag;
dtOserie(2:2:end)=duration;
time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
time=[time(2:end);time(end)];

He=exprnd(0.32,numberevents,1);
surge=ones(numberserie,1);
surge(1:2:end)=0;
Hbasesurgethrhsold=1;
surge(2:2:end)=Hbasesurgethrhsold+He;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2-Swell wave direction
angleSWELLserie=NaN*ones(numberserie,1);
    randdir=rand(numberserie,1);
    dirsign=ones(numberserie,1);dirsign(randdir<=P.Pswelldir)=-1;

    %option 1
    %rndhl=rand(numberserie,1);
    %a=find(rndhl>P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(rand(length(a),1)*45);
    %a=find(rndhl<=P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(45+rand(length(a),1)*45);

    %option 2
    angleSWELLserie=dirsign.*(rand(numberserie,1)*45/2);

%3-SeaWave direction
angleWINDserie=rand(numberserie,1)*360; %every time step a random direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%option 1
[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrierbigfarlongHOLOCENE(P);
%option 2
%load filename

%execture this command at any point, and it will save the current status of the simulation. Then you can use option 2 to initialize from here
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,'filename');


%Store value for mass balance check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumY1IN=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1IN=sum(sumY1IN(A==1));
sumY2IN=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2IN=sum(sumY2IN(A==1));
sumY3IN=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3IN=sum(sumY3IN(A==1));
FLX1=zeros(4,1);FLX2=zeros(4,1);FLX3=zeros(4,1);KBTOT=0;Y2OX=0;
FQsW_L=0;FQsW_R=0;
pondloss=0;

IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
IO.Active=Active;
fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
fIO.FQsW_R=FQsW_R;fIO.FQsW_L=FQsW_L;

zbedo=z-msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




makevideo=1;
v=VideoWriter('NAMEOFVIDEO','Motion JPEG AVI');

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')
if makevideo==1;open(v);end 
s=0;step=0;tic;
for t=1:tmax;  %iteration over the tmax time stpes
   
%Swellwave direction
angleSWELL=angleSWELLserie(t);
%SeaWave direction
angleWIND=angleWINDserie(t);%rand(1)*360; %every time step a random direction
%Lentgh of event
dtO=dtOserie(t);
%Storm surge height
Hsurge=surge(t);



if t==1;dto=0.00001;else;dto=dtO;end   
dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=P.limitdeltaz+1;maxup=P.limitmaxup+1;
        while maxdeltaz>P.limitdeltaz | maxup>P.limitmaxup
        if firstattemp==1;else;dt=dt/2*min(P.limitdeltaz/maxdeltaz,P.limitmaxup/maxup);end;firstattemp=0;
        if t<=2;dt=min(0.2*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,dx,dt,zb,IO,fIO,Hsurge,angleSWELL,angleWIND,t);
        step=step+1; %this is how many time you called the function mainevolution step
        end

    %the partial updating step was succefull! Keep going
    IO=IOtemp;
    fIO=fIOtemp;
    dti=dti+dt;%how much you have moved forward in the step
    dt=max(0,dto-dti);%the remaining time in the time step
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(t,tINT)==0;s=s+1;  
%read the variables
names = fieldnames(IO);
for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end

%read the fluxes
names = fieldnames(fIO);
for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end
    
%read the plot
names = fieldnames(PLT);
for i=1:length(names);eval([names{i} '=PLT.' names{i} ';' ]);end


z=-zb+(Yb+plyr*P.dlyr)+(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=0;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=0;
zbed=z-msl;

if time(t)>2000   
P.RSLR=1/1000/365;  %from mm/yr to m/day (the time unit is the day!)
end
% if mod(t,100)==0
% %if mod(t,250)==0
%    savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,strcat...
%        ('Barrier_slope1co20fox05W30sand_R20_t',num2str(t)));
% end


ax1 = subplot(2,2,1);%s+1 %set(IM,'alphadata',~(A==0));s+1
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-20 3],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
%colormap('jet')
caxis([-20 3]);
%colorbar('hori')
%ylim([0 60])


ax2 = subplot(2,2,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
IM=imagesc(x,y,1000*SSC2');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
colormap(ax2,'jet')
caxis([0 30]);


% ax3 = subplot(2,2,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
% %ax3 = subplot(3,1,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,flyrU2'-1);axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
% %colormap(ax3,flipud(gray))
% cptcmap('mycmap.cpt', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.01]);shading flat
% %caxis([0 1]);
% %colorbar

ax4 = subplot(2,1,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
%ax4 = subplot(3,1,3);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,Y2');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
%  colormap(ax4,'jet')
%  caxis([0 2])
 
transx_y=1;Mi=M/2; %cross section along direction x(1) or y(2)
[si,zi,Ai]=getstrat2plot(-zb,flyrU1+flyrU3,flyrb1+flyrb3,flyr1+flyr3,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
hold off;pcolorCENTER(si-dx/2/1000,zi,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
cptcmap('mycmap.cpt', 'mapping', 'direct','ncol',256); 
caxis([-1 1.01]);shading flat
%colorbar
hold on;plot(si,si*0+msl+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+msl+P.Trange/2,'--k',si,si*0+msl-P.Trange/2,'--k')
hold on;plot(si,z(:,Mi),'-k')
%VG=VEG(:,Mi);zV=z(:,Mi);hold on;
G=double(B(:,Mi)>0);G(G==0)=NaN;
plot(si,(z(:,Mi).*G),'.g')
ylim([-5 +60]);%caxis([0 1])
ylim([-5 +60]);%caxis([0 1])

% ax1 = subplot(4,1,3);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,Ux');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 0.5]);
% ax1 = subplot(4,1,4);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,Uy');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([-0.5 0.5]);


%title(strcat(num2str(time(t)),' years ',num2str(step)))
title(strcat(num2str(floor(time(t)-11600*0)),' years  (R=',num2str(P.RSLR*365*1000),'mm/yr  ',...
    '  H=',num2str(P.Ho),'m  r=',num2str(P.Trange),   'm )  #steps=',num2str(step)))



%TERM1; Qouthriver: if postive it enters
%TERM2; Qseatide: if postive it exits
%TERM3; Qseariver: if postive it exits
%TERM4 Qmouth tide. THIS IS IMPOSED ZERO BY setting D=0 at the mouth in sedtran

%SAND
QmouthRiver=FLX1(1);QseaTide=FLX1(2);QseaRiver=FLX1(3);QmouthTide=FLX1(4);
sumFLUX1=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%MUD
QmouthRiver=FLX2(1)*0;
QseaTide=FLX2(2);
QseaRiver=FLX2(3)*0;
QmouthTide=FLX2(4);
sumFLUX2=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%ORG
QmouthRiver=FLX3(1);QseaTide=FLX3(2);QseaRiver=FLX3(3);QmouthTide=FLX3(4);
sumFLUX3=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;

sumY1=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1=sum(sumY1(A==1));
sumY2=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2=sum(sumY2(A==1));
sumY3=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3=sum(sumY3(A==1));

%NOTE: Thsi is the equivalent volumetric flux, not the mass flux
checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R    [(sumY2IN-sumY2)+sumFLUX2/dx^2]-pondloss+KBTOT-Y2OX   [(sumY3IN-sumY3)+sumFLUX3/dx^2]];
if abs(checksum(1))>0.1 |  abs(checksum(2))>0.1 | abs(checksum(3))>0.1 ;checksum,pause;end

%if t==1;pause;end
if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end
pause(0.0001)
end

end

if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO



