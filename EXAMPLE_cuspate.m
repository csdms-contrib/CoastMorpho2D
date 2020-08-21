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

%Various initiliaztion for plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first line is bottom %second line is top
% %this is the good one
cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmap.cpt', cdata, ' ');

 cdata = [-1    50  255 255     0     255    255    255
          0    255 255 255     1     0    0    0];
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmapCLR.cpt', cdata, ' ');

cdata = [-1    205  165 0     0     0    0    0];
dlmwrite('mycmap1.cpt', cdata, ' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%to get always the same random numbers
rng(2)
storei=0;


%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=struct;
P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.1;%0.5;%0.2;%0.2;%0.1; % minimum water depth [m]. 0.2 is good NEEDS TO BE SMALLER THAN hwSea_lim
P.DiffSmud=1; %coefficient for tidal dispersion [-]. DO NOT CHANGE
P.DiffSsand=1;
P.DoMUD=1;%10;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)

%Sea level rise
P.RSLR=0/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%Tide
P.Ttide=12.5/24; %tidal period [day]
P.Trange=0; %mean tidal Trange [m]
P.TrangeVEG=P.Trange;%tidal Trange for vegetation [m]. Generally same of tidal range

%Storm surge
P.alpha_surge=0.25;%1;%2*0.25;%0.25;%0.25;%how much the surge contributes to tidal prsim   0.1;%

%Swell waves
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
P.Pswelldir=0.5;%0.5;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0.9; %if zero, only low angle waves
P.Ho=2;%%1.7;%3;%1.9;%1.9;%1.7;%2; %boundary swell height Hs [m]
P.Tp_swell=10;%8;%8;%6;% %boundary swell period Tp [m]
P.nrefrac=4;%either 0,1,2,3,4  Wave refraction. If zero there is no wave refraction
P.multifrequency=0;%on/off
P.wavediffraction=1;%on/off
P.Cbr=0.55;%
P.Cbed=0.038;%
P.wavefrictionCollins=0;

%Wind for sea waves
P.wind=7;%reference wind speed [m/s]

%Edge erosion
P.aw=0.3/365; %wave edge erodability m/yr/W/m2
P.maxedgeheight=2;
P.fox=0;%fraction of edge eroded material that is oxidized.

%Wind waves and swell numerics
P.hwSwelltransport_lim=0.5;%0.2;%1m  %after April 16th 2019 this limite the along-shore transport (from both alongwave and radiation stress)!!!
P.hwSwell_lim=0.10001;%2;%0.2 %limiter water depth for swell
P.hwSea_lim=0.2;%0.10001;%.2;%0.5;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE" ; NEEDS TO BE LARGER THAN KO!!!!

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=0*10/1000; %Sea boundary SSC for mud [g/l]
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
P.ws1=0.02/1;%sand with D50=500um  0.05 %m/s
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);

%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%
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
P.fMFswell=1;%0.06; %use 1 if you use the equvilanet wave height
P.fMFsea=1; %use 1 if you use the equvilanet wind speed
P.fMFriver=10/365;

%Vegetation parameters
P.dBlo=-(0.237*P.TrangeVEG-0.092);
P.dBup=P.TrangeVEG/2;%-0.2;
P.Cv=0.1;%Manning for vegetated ara
P.wsB=1/1000;%Mud Settling velocity for vegetated ara
P.taucrVEG=0.5;%Critical sheak stress for vegetated areas

%Organic accretion by vegetation
P.AccreteOrganic=0;
P.Korg=0/1000/365;%5/1000/365;%P.Korg=org/365;%5/1000/365;

%ON/OFF processes
P.VEGETATION=0;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.depthlimiterflow_withVEG=0;
P.computemud=0;
P.computesand=1;
P.computeSwellwave=1;
P.computeSeaWaves=0;
P.computeEdgeErosionSwell=0;
P.computeEdgeErosionSea=0;
P.compute_currentbankerosion=0;
P.computetide=0;
P.computeriver=0;

%Correction for second-order river dynamics
P.riverwaterlevel=0;
P.rivermomemntumcorrection=0;

%Various boundary condtions
P.periodic=1;
P.imposeseaboundarydepthmorphoALL=0; %to use when a channel mouth is at a boundary

%Ebb-flood momentum correction
P.ebbfloodcorrection=1;
P.residualcurrents=0;

%Curvature flow modifications
P.curvaturecorrection=0;
%P.alphacurvaturesmooth=100;
%P.curvfactor=40;

%Pond dynamics
P.calculateponddynamics=0;
P.parameterforpond1=1;
P.parameterforpond2=1;
P.parameterforpond3=1;

%Stratigraphy
P.evolvestratigraphy=0;
P.VEGstratigraphy=0;%if 0 then you put the organic into the mud. If 1 then you calculate the organic as a sediment per se (advection, divergence,etc)
P.VEGonsand=0;

%Stratigraphy parameters
P.conces=10;%how much to extra erode, a parameter
P.nlyr=20; %max number of layers
P.dlyr=1; %thickenss of layers
P.tlyrU=3; %max depth to add layer %must be larger than dlyr
P.tlyrD=0.5; %min depth merge layers %mus be larger than dlyr
P.tcko=10;%tickness of bed layer
P.levo=15;%intial level occupied
P.YUi=1000;%initial thickess of active layer
P.initialfU=1;%initial composition of the active layer
P.initialf=1;%initial composition of all the layers


P.reducefractionsediment=1;%this should be 1 unless you to strange stuff. ADDED JUNE 2019@@@@@@@@@@@

%Global numerical limiters
% limitdeltaz=10;
% limitmaxup=5;
limitdeltaz=10/2;
limitmaxup=5/2;
P.limitdeltaz=limitdeltaz;
P.limitmaxup=limitmaxup;


%Time parameters
tmax=20000;%4000;%3000;%13385;%15000;%200;%110;%1250;% number of time steps
tINT=1;%how many time steps you want to do the plot (if 1 you plot every time step). Does not affect the computation
%dtO=2*365;%the reference time step [days]. The reference unit of time is days!!!
%time=[0:tmax-1]*dtO/365; %converted to years, just to plot. Does not affect the computation




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time series
numberserie=20000;%if you change this you will change the actual values in the time series, rememebr!
numberevents=numberserie/2;
% lag=exprnd(1.7*365,numberevents,1);lag(lag<1)=1;
% duration=exprnd(0.01*365,numberevents,1);lag(lag<0.01)=0.01;
% dtOserie=ones(numberevents*2,1);
% dtOserie(1:2:end)=lag;
% dtOserie(2:2:end)=duration;
% time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
% time=[time(2:end);time(end)];
% 
% He=exprnd(0.5,numberevents,1);
% % Hoseries=ones(numberevents*2,1);
% % Hoseries(1:2:end)=2;%1;
% % Hoseries(2:2:end)=4+He;
% 
% surge=ones(numberevents*2,1);
% surge(1:2:end)=0;
% %surge(2:2:end)=0.5+He*0.6;%+(-0.5+rand(numberevents,1));
% surge(2:2:end)=0.5+He*1;%+(-0.5+rand(numberevents,1)); 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



lag=exprnd(1.8*365,numberevents,1);%lag(lag<1)=1;
duration=5*ones(numberevents,1);%exprnd(2,numberevents,1); %duration=exprnd(0.26,numberevents,1);lag(lag<0.01)=0.01; %duration=exprnd(0.0003*365,numberevents,1);lag(lag<0.01)=0.01;
dtOserie=ones(numberserie,1);
dtOserie(1:2:end)=lag;
dtOserie(2:2:end)=duration;
dtOserie(:)=0.05*365;%%%%%%%%%%%%%
%dtOserie(:)=5*365;%%%%%%%%%%%%%
time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
time=[time(2:end);time(end)];

He=exprnd(0.16,numberevents,1);
% Hoseries=ones(numberevents*2,1);
% Hoseries(1:2:end)=2;%1;
% Hoseries(2:2:end)=4+He;

surge=ones(numberserie,1);
surge(1:2:end)=0;
%surge(2:2:end)=0.5+He*0.6;%+(-0.5+rand(numberevents,1));
%Hsuregwaverunup=0.5;
Hbasesurgethrhsold=1;
surge(2:2:end)=Hbasesurgethrhsold+He;%+(-0.5+rand(numberevents,1)); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%Swell wave direction
angleSWELLserie=NaN*ones(numberserie,1);
    randdir=rand(numberserie,1);
    dirsign=ones(numberserie,1);dirsign(randdir<=P.Pswelldir)=-1;
rndhl=rand(numberserie,1);
    a=find(rndhl>P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(rand(length(a),1)*45);
    a=find(rndhl<=P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(45+rand(length(a),1)*45);
%angleSWELLserie=dirsign.*(rand(numberserie,1)*45/2);

%SeaWave direction
angleWINDserie=rand(numberserie,1)*360; %every time step a random direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrierBARRIER(P);
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrierbigfarlong(P);
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrierbigfar(P);
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrierbigfarlongHOLOCENE(P);

[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,S]=initializegeometry_3sedimentscuspate(P);


%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell]=initializegeometry_3sedimentsbarrier_restoration(P);

%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S,'SEGA_R20R1_H17T8long_newsurge_aeolian');
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S,'SEGAR1steadystate');
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,S,'SEGAturbo_R10');
%SEGA_R1_H17T8long_newsurge_aeolian_dx100
%load HOLOCENEtemp11309
%load HOLOCENEall
%load wavetranslim1_anglespread45_long_R3
%load G
%load SEGA_R1_H17T8long_newsurge_aeolian_dx100

%Y1(197:200,80:110)=Y1(197:200,80:110)+3;
%Y1(215+4:218+4,100:130)=Y1(215+4:218+4,100:130)-5;
%Y1(230:233,60:100)=Y1(230:233,60:100)+5;
%Y1(230+4:233+4,60:100)=Y1(230+4:233+4,60:100)-5;
%load SEGAturbo_R10

%load Gang22R1
%load basin_noangleokbisBBBfaclong1bis
%load Basin
% Initck=2;
% zb=zb-Y1+Initck;
% Y1=Y1*0+Initck;
%load BhlimVEG___MUD
%zb(end-5:end,:)=zb(end-5:end,:)+50;



%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,'cuspatedx100');
%load cuspatedx100
%Active(:)=1;

%barrier restoration
% for i=80:120
%     a=find(zbedo(:,i)<0.5);a=a(end);
%     Y1(a-1:a+6,i)=Y1(a-1:a+6,i)+4;
%     %Y1(a-1+8:a+6+8,i)=Y1(a-1+8:a+6+8,i)-4;
%     Y1(240-1:240+6,i)=Y1(240-1:240+6,i)-4;
% end
%zbedo=z-1000;
% for i=180:230
%     a=find(zbedo(:,i)>-4);a=a(end);
%     a=a-10;
%     Y1(a-1:a+6*2,i)=Y1(a-1:a+6*2,i)+4;
%     %Y1(a-1+14:a+6*2+14,i)=Y1(a-1+14:a+6*2+14,i)-4;
%     Y1(a-1+14+80:a+6*2+14+80,i)=Y1(a+80-1+14:a+6*2+14+80,i)-4;
%     %Y1(240*2-1:240*2+6*2,i)=Y1(240*2-1:240*2+6*2,i)-4;
% end

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



makevideo=0;
%v=VideoWriter('FRIDAYdx200_R20_downslopewave2_long_xxx','Motion JPEG AVI');
v=VideoWriter('Cuspate_U07_dwslopex1VEROMAX_Qw02_COS_long39_lonthr1m_slopeMAX','Motion JPEG AVI');

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')
if makevideo==1;open(v);end 
s=0;step=0;tic;
for t=1:10000%2500%tmax;%11309:tmax;  %iteration over the tmax time stpes
    
%if t>4120;P.RSLR=1/1000/365;end
%if t>1157;P.RSLR=3/1000/365;end
%if t>200;P.Ho=0;end

%Swellwave direction
angleSWELL=angleSWELLserie(t)
%SeaWave direction
angleWIND=angleWINDserie(t);%rand(1)*360; %every time step a random direction
%Lentgh of event
dtO=dtOserie(t);
%Storm surge height
Hsurge=0;%surge(t);

% if Hsurge>0;P.Ho=2*1.7;%3;%1.9;%1.9;%1.7;%2; %boundary swell height Hs [m]
% else   
% P.Ho=1.7;%3;%1.9;%1.9;%1.7;%2; %boundary swell height Hs [m]
% end

[dtO P.Ho Hsurge]
if t==1;dto=0.00001;else;dto=dtO;end   
dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;
        if t<=2;dt=min(0.2*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,dx,dt,zb,IO,fIO,Hsurge,angleSWELL,angleWIND,t);
        step=step+1; %this is how many time you called the function mainevolution step
        end

    %the partial updating step was succefull! Keep going
    IO=IOtemp;
    fIO=fIOtemp;
    dti=dti+dt;%how much you moved forward
    dt=min(dt*2,max(0,dto-dti));%the remaining time in the time step
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



if mod(t-1,5)==0
    storei=storei+1;t
    STORE(:,:,storei)=zbed;
end

ax1 = subplot(1,2,1);%s+1 %set(IM,'alphadata',~(A==0));s+1
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-10 10],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
%colormap('jet')
caxis([-10 10-0.2]);
%colorbar('hori')
%ylim([0 60])


ax2 = subplot(1,2,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(x,y,Hs');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
colormap(ax2,'jet')
caxis([0 P.Ho]);

%title(strcat(num2str(time(t)),' years ',num2str(step)))
title(strcat(num2str(floor(time(t)-11600*0)),' years  (RSLR=',num2str(P.RSLR*365*1000),'mm/yr)   #steps=',num2str(step)))



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
checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R    [(sumY2IN-sumY2)+sumFLUX2/dx^2]-pondloss+KBTOT-Y2OX    [(sumY3IN-sumY3)+sumFLUX3/dx^2]];
if abs(checksum(1))>0.1 |  abs(checksum(2))>0.1 | abs(checksum(3))>0.1 ;checksum,pause;end


if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end
pause(0.1)
end

end

if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO



