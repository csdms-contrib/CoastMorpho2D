clear; close all;clc


rng(2)%to get always the same random numbers

figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')
for ccc=1;

cd(fileparts(matlab.desktop.editor.getActiveFilename));
[FILEPATH]=fileparts(pwd);
%[FILEPATH]=fileparts(FILEPATH);
addpath(FILEPATH);


%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ON/OFF processes
P.VEGETATION=0;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.computemud=0;
P.computesand=1;
P.computeSwellWave=0;
P.computeSeaWave=0;
P.computeEdgeErosionSwell=0;
P.computeEdgeErosionSea=0;
P.computetidalcurrent=0;
P.computeRiver=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%General Hydrodynamics - numerics
P.g=9.81; %gravity [m/s2]
P.kro=0.01;%0.01;%0.1;%0.2;%1;%0.2;%0.02;%02;%0.5;%0.2;%0.2;%0.1; % minimum water depth [m]. 0.2 is good NEEDS TO BE SMALLER THAN hwSea_lim
P.hlimC=0.01;%0.5;%limiter to calcualte SAND flow erosion 
P.ftidepeak=5;%25;%2/(Ttide/2*24);%1 hour every half tide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%General Hydrodynamics - site specific
%Sea level rise
P.RSLR=0/1000/365;  %from mm/yr to m/day (the time unit is the day!)
P.COMP=0/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%make the cells next to A==0 extra firction
P.solidboundaryextrafriction=0;

%eddy diffusivity + no slip
P.ahBNK=0.5;
P.accumulateonedges=0;

%Manning coeffinent unvegeated (same for sand and mud)
P.CbT=0.02;
P.FcrUT=0.3;
P.FcrUTveg=0.3;%0.1;
P.UTref=0.5;

%Manning for sediment transport
P.CbS_MUD=0.02;
P.CbSveg_MUD=0.03;
P.CbS_SAND=0.02;%0.02;%0.015;%0.022;%0.015;%0.022;%15;
P.CbSveg_SAND=0.03;

%Variable tide
P.tidalrangeattenuation=0;
% P.Trange90_o=0.51;%these remain constan as boundary conditions
% P.Ttide90=1.02;
% P.MSL90=0.15;

%Vegeation limits
P.trackmarshage=0;
P.dBlo=-0.1;
%P.Dynamicvegationrange=1;
P.Dynamicvegationrange=0;
%P.dBup=(0.37+0.1);%(2.06+0.1)/2;%-0.2;%P.TrangeVEG=(2.06+0.1)*2;%P.Trange;%P.Trange;%tidal Trange for vegetation [m]. Generally same of tidal range


% %Swell waves parameters
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
% P.nrefrac=0;%either 0,1,2,3,4  Wave refraction. If zero there is no wave refraction
% P.multifrequency=0;%on/off
% P.wavediffraction=1;%on/off
% P.Cbr=0.73;%73;
% P.Cbed=0.038;%NaN;%wave bed friction if you use Jonswap (0.067 or 0.038)
% P.wavefrictionCollins=0;

%Wind waves and swell numerics
P.hwSwelltransport_lim=1;
P.hwSwell_lim=0.2; %limiter water depth for swell: swells are imposed zero below this limit
P.hwSea_lim=0.2;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE"%NEEDS TO BE LARGER THAN KO!!!!
P.extraHseaimposed=0;

%Shallow tidal flow
P.calculateshallowflow=0;
P.kroS=0.01;
P.nSHALLOW=10;
P.maxQshallow=1;
P.FcrUS=0.3;

%Subaerial: aeolina and runup
P.Procesaeolian=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%General Sediment
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.DiffSsand=0.5;%0.2; %coefficient for tidal dispersion [-]. 1 DO NOT CHANGE
P.DiffSmud=1;%0.5;%0.2; %coefficient for tidal dispersion [-]. 0.1 DO NOT CHANGE
P.DoMUD=10;%10;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)
P.DoSAND=0.1;%1;%0.1;%0.1;

%Edge erosion
P.aw=0.1/365; %wave edge erodability m/yr/W/m2
P.maxedgeheight=2;
P.fox=0.25;%fraction of edge eroded material that is oxidized.

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=50/1000;%40/1000; %Sea boundary SSC for mud [g/l]
P.co3=0/1000; %Sea boundary SSC for mud [g/l]
P.SeaSSCbelowLIMIT=1;
P.ZlevcoLIMIT=-1;%[m]

%Sand
P.d50_1=0.25/1000;% %sand grain size [m]
P.ws1=0.02*1*1;%       (sand with D50=500um will have 0.05 %m/s)
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);

%Sand Parameters: Downslope paramters for sand and mud (proportional to the sediment transport Qs!!!!)
P.alphaSAND=NaN; %coefficient for bedload downslope of sand. Calibrated with JMSE 2018, do not change
P.alphaSANDriver=1;%2;%0.1;%0.5;%3;%2;%5;%0.5; %coefficient for bedload downslope of sand. Calibrated with JMSE 2018, do not change
P.hlimCdwnSAND=0;%limit to apply to downslope of sand (and maybe also mud) Calibrated with JMSE 2018, do not change
P.downslopeSANDseawaves=NaN;%2;%10;%0;%5*10;% multiplication for sea-waves dowbsloep transport (^5 as for the swells) for SAND
P.reduceSANDbankVEG=0.1;

% if ccc==1;P.alphaSANDriver=0.1;
% elseif ccc==2;P.alphaSANDriver=1;
% elseif ccc==3;P.alphaSANDriver=2;
% elseif ccc==4;P.alphaSANDriver=5;
% elseif ccc==5;P;P.alphaSANDriver=10;
% end
    
    P.alphaSANDriver=5;
    
%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%0.2/1000;% m/s
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);%P.por2=0.7
%P.por2=0.8;P.rbulk2=P.rhos*(1-P.por2);

%Mud parameters
P.me=0.15*10^-4*24*3600;  %per day!!!
P.taucr=0.2;
P.tcrgradeint=0;% Pa/m
P.leveltauincrease=NaN;%P.TrangeVEG/2;%1;
P.crMARSH=0.1/365;%creep coefficient vegetated
P.crbank=0.1/365;%0.5/365;%creep coefficient vegetated
P.crMUD=3/365;%creep coeffcinet
P.alphaMUD=0.2; %coefficient for bedload downslope of mud. added April 2019. Similar to P.alphaSAND. Changed to 0.25 from 0.5 because initially used bulk1 instead of bulk2
P.facQsbank=0.001*365/365;
P.facHwave=0.5;
P.hlimCdwnMUD=1;
P.crSAND=0/365;


%%bank erosion parameters
P.hDlo=0.01;%NOT RELEVANT IF pDo=pDm
P.hDhi=0.1;%NOT RELEVANT IF pDo=pDm

P.pDo=1;%0.5;
P.pDm=1;%0.5;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vegetation parameters
P.Cv=0.05;%Manning for vegetated ara
P.wsB=0.2/1000;%Mud Settling velocity for vegetated ara
P.taucrVEG=0.2;%Critical sheak stress for vegetated areas

%Dynamic vegeation growth
P.Dynamicvegeationgrowth=0;
P.bioseed=0.01;
P.dBseed=P.dBlo;
P.plantexpansionrate=1;
P.abiogrow=1/(1*365); %unit of 1 over time. The larger the denominator, the slower it grows

%Organic accretion by vegetation
P.AccreteOrganic=1;
P.Korg=5/1000/365;%8/1000/365;%P.Korg=org/365;%5/1000/365;
P.variableORGaccretion=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pond dynamics
P.calculateponddynamics=0;
P.variablePONDdynamic=0;
    P.Epondform=1*10^-4/100;
    P.zpondcr=-0.2;%P.Trange/4;%base of new pond formation with respect to MSL
    P.minponddepth=0.1;%1; %minimum depth to define a pond after you identified the "lakes"
    P.maxdpond=max(0.2,max(P.minponddepth*1.2,0.3));%0.5;%0.5;%maximum depth scour of a new pond
    %%
    P.zntwrk=0;%(P.Trange/2)*0.5;%(P.Trange/2)*0.2;%P.Trange/2*0.9;%P.Trange/2-0.3;%0.3; %depth above msl that defines the channel network.  the smaller the harder to drain!
    P.distdr=NaN;%Clealry if nan is nto used4; %m, distance of extra drainage
    %%
    P.aPEXP=0.05;%0.015*10;%isolated pond expansion rate m/yr
    P.ponddeeprate=0.003;%m/yr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%River inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.fMFriver=1;%100/365;%Correction for proceeses duration (to scale waves and tidal transport, because waves do not occur all the time!)
P.CbR=0.02;
P.Uref=1;
P.minUfric=0.2;%0.2;%0.5;
P.nonlinearfriction=1;P.URfrictiono=NaN; %this is the reference velcoity used if linear friction
P.imposenocrossflux=1;%This imposes a uniform velocity at the river mouth. 
P.imposeNOerosiondepostionatmouthSAND=0;
P.imposeNOerosiondepostionatmouthMUD=0;
P.imposeSANDCeqartivermouth=0;
P.AcrreteRivermouthwithRSLR=1;
%Impose the sediment discharge input at the river mouth
%INPUT a: A==10 -Davis Pond
%hmouth_a=10; %water depth [m]%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
%Qmouth_a=hmouth_a*2.1; %river discharge per unit of cell [m2/s]
hmouth_a=2.6; %water depth [m]%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
Qmouth_a=4*1;%hmouth_a*1.54;%1.32; %river discharge per unit of cell [m2/s]
co2mouth_a=NaN;%100/1000;%500/1000; %SSC of mud at the river [g/l]
directionQ_a=1; %1:up -1:down 2:left  -2:right

% hmouth_b=5; %water depth [m]%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
% Qmouth_b=5*0.5; %river discharge per unit of cell [m2/s]
% co2mouth_b=100/1000;%500/1000; %SSC of mud at the river [g/l]
% directionQ_b=-2; %1:up -1:down 2:left  -2:right

P.Qmouth=[Qmouth_a];% Qmouth_b];
P.hmouth=[hmouth_a];% hmouth_b];
P.co2mouth=[co2mouth_a];%  co2mouth_b];
P.co1mouth=[100/1000]*1;% 0];
P.directionQ=[directionQ_a];% directionQ_b];


%Correction for second-order river dynamics
P.FcrUR=NaN;%0.3;
P.FcrURflow=1;
P.riverwaterlevel=1;
P.nitero=10;
P.initialhforriverflow=1;%[m]%
P.usepreviousFLOW=1;P.hpRIVo=[];


P.errorWL=0.2;P.errorUUfac=P.errorWL/0.2;P.checkerrorh=0.2;P.checkerrorU=0.2;P.factorRiverIter=0.05;P.hscaleforiteration=5;
P.hITERminBELOWbed=-0.5;
P.hITERupMAX=0.1;
P.NITMAX=100;

P.correctmomentumriver=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Various boundary condtions
P.periodic=0;
P.wavetransportzerolateralgradient=0;%only needed for waves, if period=0
P.imposeseaboundarydepthmorphoALL=0; %to use when a channel mouth is at a boundary


%Tidal non-linear friction
P.tidalnonlinearflow=1;

%Ebb-flood momentum correction
P.correctmomentumtideebbflood=0;
P.residualcurrents=0;

%CorrectMomentum parameters(for both tide and river)
P.DDUlimit=0.99;%0;%.5;
P.facmomentumRIVER=0.2;
P.facmomentumTIDE=0.05;

%Curvature flow modifications
P.curvaturecorrection=0;
P.flowbankerosion=1;%This enables the process related to a_bankerosion
P.a_diffusecur=0.01;%smooth the curvature in the transverse direction. Not very sensitive to this. Just keep is high (>0.1)
P.advectflow=30;%200;%300;%250;%3000;%Modify the flow (not change. if change, only 10-20%). Already include dx factor in mainevolutionstep.
P.a_bankerosion=0.1/365;%active bank erosion - can be changes (values 1-100)
P.a_bankerosionUNVEG=0.1/365;

%Connect channel corners to make channels more continous
P.connectchannelcorners=0;

%Caluclate salinity
P.calculatesalinity=0;
P.SALTocean=25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stratigraphy
P.evolvestratigraphy=0;
P.VEGstratigraphy=0;%if 0 then you put the organic into the mud. If 1 then you calculate the organic as a sediment per se (advection, divergence,etc)
P.VEGonsand=0;
%Stratigraphy parameters
P.conces=10;%how much to extra erode, a parameter
P.nlyr=20; %max number of layers
P.dlyr=0.5; %thickenss of layers
P.tlyrU=0.8; %max depth to add layer %must be larger than dlyr
P.tlyrD=0.2; %min depth merge layers %mus be larger than dlyr
P.tcko=10;%tickness of bed layer
P.levo=15;%intial level occupied
P.YUi=1000;%initial thickess of active layer
P.initialfU=1;%1;%0.5;%initial composition of the active layer
P.initialf=1;%1;%0.8;%initial composition of all the layers

P.reducefractionsediment=1;%this should be 1 unless you to strange stuff. ADDED JUNE 2019@@@@@@@@@@@
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Global numerical limiters
limitdeltaz=2;%10;%/4;
limitmaxup=0.1;%0.5*4;%5/4;%0.5;%

%Time parameters
tmax=2*501;%1500*1;%10000;%2000;%400;%2*199;%1000;%149;%1000;%149;%2000;%3000;%50;%42/2-1;%2000;%1250;% number of time steps
tINT=10;%how many time steps you want to do the plot (if 1 you plot every time step). Does not affect the computation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%time series%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numberserie=15000;%2000;%2000;%10000;%if you change this you will change the actual values in the time series, rememebr!
numberevents=numberserie/2;
lag=exprnd(1.7*365,numberevents,1);lag(lag<1)=1;
duration=exprnd(0.01*365,numberevents,1);lag(lag<0.01)=0.01;
dtOserie=ones(numberevents*2,1);
dtOserie(1:2:end)=lag;
dtOserie(2:2:end)=duration;


dtOserie=dtOserie*0+0.01*365/1*2;% /2;
dtOserie(1:20)=0.01*365*1/5;

time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
time=[time(2:end);time(end)];

He=exprnd(0.5,numberevents,1);

surge=ones(numberevents*2,1);
surge(1:2:end)=0;
surge(2:2:end)=0.5+He*0.6;%+(-0.5+rand(numberevents,1));
surge=surge*0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Swell wave direction
P.Pswelldir=0.5;%0.5;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0.0; %if zero, only low angle waves
P.Ho=3;%1.7;%2; %boundary swell height Hs [m]
P.Tp_swell=8;%8;%6;% %boundary swell period Tp [m]
t_SWELLdir=NaN*ones(numberserie,1);
    randdir=rand(numberserie,1);
    dirsign=ones(numberserie,1);dirsign(randdir<=P.Pswelldir)=-1;
t_SWELLdir=dirsign.*(rand(numberserie,1)*45/2);

%SeaWave direction
angleWINDserie=rand(numberserie,1)*360; %every time step a random direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,B]=initializegeometry_example3DeltaLARGE(P);
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,B,x,y,msl,SPCLcell,'FILE');

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
IO.zb=zb;
IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
IO.Active=Active;
IO.B=B;
IO.AGE=A*0;
fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
fIO.FQsW_R=FQsW_R;fIO.FQsW_L=FQsW_L;

z=-zb+(Yb+plyr*P.dlyr)+(Y1+Y2+Y3);
zbedo=z-msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%track sediment flux (optional)
P.tracksedimentfluxes=0;P.XX=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




makevideo=1;
v=VideoWriter('VIDEO2_wl02','Motion JPEG AVI');

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if makevideo==1;open(v);end 
s=0;step=0;tic;
for t=1:tmax  %iteration over the tmax time stpes
   

dtO=dtOserie(t);

%Not used but needed to be defined
P.Trange_o=0;
P.tempdeltaMSL=0;
P.WINDspeed=0;
P.WINDdir=0;
P.angleSWELL=0;
P.Ttide=NaN;


if t==1;dto=0.00001;else;dto=dtO;end   


dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;
        if t<=2;dt=min(0.5*365,dt);end
        %[IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstepRIVERCASE(A,AW,SPCLcell,P,x,dx,dt,IO,fIO);
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,x,dx,dt,IO,fIO);
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
if t>0%==149;%>0;%==249;
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


ax1=bigsubplot(1,1,1,1,0.147,0.005);%,[hgap],[vgap])
IM=imagesc(y,x,zbed );set(IM,'alphadata',~(A==0));set(gca,'YDir','reverse');%colormap('jet');
cmp=demcmap([-5 2],256); %-3 1 P.Trange/2
colormap(ax1,cmp);caxis([-5 2]);axis equal
title(strcat(num2str(time(t)),' years ',num2str(step)))
text(1,0.1,strcat(num2str(time(t)),' years ',num2str(step)))




%SAND
QmouthRiver=FLX1(1);QseaTide=FLX1(2);QseaRiver=FLX1(3);QmouthTide=FLX1(4);
sumFLUX1=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%MUD
QmouthRiver=FLX2(1);QseaTide=FLX2(2);QseaRiver=FLX2(3);QmouthTide=FLX2(4);
sumFLUX2=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%ORG
QmouthRiver=FLX3(1);QseaTide=FLX3(2);QseaRiver=FLX3(3);QmouthTide=FLX3(4);
sumFLUX3=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;

sumY1=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1=sum(sumY1(A==1));
sumY2=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2=sum(sumY2(A==1));
sumY3=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3=sum(sumY3(A==1));


checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R    [(sumY2IN-sumY2)+sumFLUX2/dx^2]-pondloss+KBTOT-Y2OX    [(sumY3IN-sumY3)+sumFLUX3/dx^2]];
if abs(checksum(1))>10 |  abs(checksum(2))>10 | abs(checksum(3))>10 ;checksum,pause;end

if (makevideo==1 & mod(t,2)==0) ;V=getframe_nosteal_focus(1,[1600 900]);writeVideo(v,V);end
pause(0.01)
end

end


end; %end of panel plotting
if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO

end
toc