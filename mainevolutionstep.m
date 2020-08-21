function [IO,fIO,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,PARAMS,dx,dt,zb,IO,fIO,Hsurge,angleSWELL,angleWIND,t)


%Load the parameters and variables
names = fieldnames(IO);
for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end

names = fieldnames(fIO);
for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end

names = fieldnames(PARAMS);
for i=1:length(names);eval([names{i} '=PARAMS.' names{i} ';' ]);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(A);
optionBC=2; %impose water depth (for sea boundary)

%note: zb has the opposite convention, is positive going down! That is why
%you have to deal with -zb. Just think a zbedrock=-zb;

Y=Y1+Y2+Y3;%total thickness of active later
zs=-zb+(Yb+plyr.*(dlyr)); %the heigth of the ground just below Y
%%%zs=zberock+(Yb+plyr.*(dlyr)); %the heigth of the ground just below Y
z=zs+Y;%the absoulte bed elevation, positive is the bed is high ground, negative if it is low ground

zoriginal=z;

%Volumetric fraction
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
if reducefractionsediment==1;
fracY1=max(0,Y1)./Ytot;fracY1(Ytot==0)=0;%volumetric fraction of sand in active layer
fracY2=max(0,Y2)./Ytot;fracY2(Ytot==0)=0;%volumetric fraction of mud in active layer
fracY3=max(0,Y3)./Ytot;fracY3(Ytot==0)=0;%volumetric fraction of organic in active layer
else
fracY1=A*0+1;fracY2=A*0+1;fracY3=A*0+1;
end
    


%Hydrodynamic%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sea level
msl=msl+RSLR*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the upland cells active
Active((z-msl)<Trange/2)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Pond dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need to run it first because it affects the vegetation through S !!!!!
if calculateponddynamics==1    
zsill=Trange/2;  
[deltaY2,pondloss]=pondformation(A,dx,dt,Epondform,z-msl,zpondcr,maxdpond,zsill,pondloss);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
[S,AC,DIF]=findisolatedponds(z-msl,A,N,M,dx,zntwrk,zsill,distdr,minponddepth);%AC is only used for plotting, not in the computation   
%DIF is the amount of water impunded at low water. It is the remainng water depth in the pond!
[S,deltaY2,pondloss]=isolatedpondexpansion(z-msl,S,A,N,M,dx,dt,zpondcr,maxdpond,aPEXP,pondloss);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
[deltaY2,pondloss]=isolatedponddeepening(z-msl,S,ponddeeprate,dt,pondloss);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
else;AC=A*0;S=A*0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Accrete in the river mouth to maintain same water depth
Y1(A==10)=Y1(A==10)+RSLR*dt; %sediment on the bed of the river mouth!!
z(A==10)=zs(A==10)+(Y1(A==10)-Y2(A==10));  %accrete the mouth ourlet with sand

%Water depth
if computeriver==1 & riverwaterlevel==1
hpRIV=max(0,z+5);% %at the beginning, just a small water level from the river to avoid extra slopes
else;hpRIV=0;end
[h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);

%Elevation with respect to MSL
Zlev=z-msl;

%Vegetation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if VEGETATION==1;
VEG=Zlev>dBlo;
B=4*(Zlev-dBup).*(dBlo-Zlev)/(dBlo-dBup)^2;B(Zlev>dBup)=0;B(Zlev<dBlo)=0;  %B=(dBup-lev)/(dBup-dBlo);B(lev>dBup)=0;B(lev<dBlo)=0;
else;VEG=A*0;B=A*0;end
B=B.*(S==0);%no biomass where there are ponds
VEG=VEG.*(S==0);%no vegeation where there are ponds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Manning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MANN=0*A+Cb;
MANN(VEG==1)=Cv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%River flow
if computeriver==1;
    Umouth=NaN; %not used anymore
    Uo=A*0+1;%first attempt of velocity
    [UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,max(0,h),dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
        if riverwaterlevel==1
            niter=5-floor(rand(1)+0.5);
            for i=1:niter
            hpRIV=(hpRIV+hriver)/2;
            [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV); 
            Uo=(Uo+UR)/2;
            [UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,h,dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
            end
        else;hpRIV=A*0;end
        hpRIV=hriver;
else;UR=0*A;URy=0*A;URx=0*A;end


%CORRECTION RIVER FOR MOMENTUM ADVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeriver==1 & rivermomemntumcorrection==1
[Vx,Vy]=correctvelocity4MOMENTUMebbflood(N,M,periodic,h,A,MANN,dx,URx,URy);
    URX=URx;a=find(URX.*Vx>0);URX(a)=URX(a)+Vx(a);
    URY=URy;a=find(URY.*Vy>0);URY(a)=URY(a)+Vy(a);
    URx=URX;URy=URY;
    UR=sqrt(URX.^2+URY.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if riverwaterlevel==0
h(A==10)=msl-z(A==10);% the river is not influenced by tide or Hsurge (these two!)
end

%Pre-calculate sediment input at river mouth
if computeriver==1;
[E,Ceq]=totalsedimenterosionSANDsine(h(A==10),hlimC,0,ss,d50_1,ws1,1,URx(A==10),fMFriver,kro,MANN,VEG,U*0);          
[QsmouthSAND]=Ceq.*h(A==10).*URx(A==10);
QsmouthMUD=URx(A==10)*co2mouth*fMFriver.*h(A==10);
QsmouthSAND=sum(QsmouthSAND);
QsmouthMUD=sum(QsmouthMUD);
else;QsmouthSAND=0;QsmouthMUD=0;end



%Tide&Surge flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetide==1;

if calculateponddynamics==1;
%DIF is the impounded water depth, need to subtract to the imput discharge from the ponded area (THE WATER REMAINS THERE!!!)
pondleakage=0.2;
DIF=max(DIF,0); %you cannot impound a negative water depth!!! This happens because of the trick used to swap the cell during the pond floodin
DIFo=DIF;
dtideI=Trange/2-(z-msl)-DIFo;
dtide(S==1)=max(0, Trange/2-(z(S==1)-msl)-max(0,DIFo(S==1)-pondleakage));
fTide(S==1)=0.01;
end
DHeff=min(Trange,dtide)+alpha_surge*min(Hsurge,dsurge);%the water level excusrion for the total water prismmin

[U,Uy,Ux]=tidalFLOW(A,MANN,h,ho,dHW,dx,DHeff,Ttide,periodic,kro);
U(A==10)=0;Ux(A==10)=0;Uy(A==10)=0;
else;U=0*A;Uy=0*A;Ux=0*A;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%CORRECTION FOR MOMENTUM ADVECTION$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ebbfloodcorrection==1
[Vx,Vy]=correctvelocity4MOMENTUMebbflood(N,M,periodic,h,A,MANN,dx,Ux,Uy);
    %EBB
    VebbX=Ux;a=find(VebbX.*Vx>0);VebbX(a)=VebbX(a)+Vx(a);
    VebbY=Uy;a=find(VebbY.*Vy>0);VebbY(a)=VebbY(a)+Vy(a);
    Vebb=sqrt(VebbX.^2+VebbY.^2);
    %FLOOD
    VfloodX=-Ux;a=find(VfloodX.*Vx>0);VfloodX(a)=VfloodX(a)+Vx(a);
    VfloodY=-Uy;a=find(VfloodY.*Vy>0);VfloodY(a)=VfloodY(a)+Vy(a);
    Vflood=sqrt(VfloodX.^2+VfloodY.^2);
        
    %residual currents
    if residualcurrents==1;
    UresX=-(VebbX+VfloodX);UresY=-(VebbY+VfloodY);
    UresX(A==2 | A==0)=0;UresY(A==2 | A==0)=0;
    UresX(1:2,:)=0;UresX(end-1:end,:)=0;%to conserve mass at the open boundary. Need to get velocity zero next to A==2 (upwind!)
    else;UresX=0;UresY=0;
    end
else;Vebb=U;Vflood=U;UresX=0;UresY=0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Swell waves%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSwellwave==1;        
    hwave=ho;%just to isolate this water depth and not mess up. It can theoretically be different than ho    
   
    %THIS IS THE OLD SINGLE FREQUENCY AND SINGLE DIRECTION
         % kwave=wavek(1/Tp_swell,hwave);
         % [Hs]=SwellWaves(A,AW,Ho,N,M,hwave,Tp_swell,kwave,dx,periodic,angleSWELL,gridDIR);  
    if multifrequency==1; %Case multiple frequency
        [Tperiodi Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,[1 1.5 2 2.5]);
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap,nrefrac,wavediffraction);
    elseif multifrequency==0 %Case single frequency  
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirection(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tp_swell,dx,periodic,angleSWELL,gridDIR,nrefrac,wavediffraction);
    end
    waveANGLE(isnan(waveANGLE))=0;

          if computesand==1;
          [QsWslope,QsWon]=WaveSedimentTransport(Hs,hwave,kwave,rhos,N,M,wavePERIOD,dx,ss,ws1,hwSwell_lim,fTide);
          QsWon(hwave<=hwSwell_lim)=0;QsWslope(hwave<=hwSwell_lim)=0;deltaPW(hwave<=hwSwell_lim)=0;
          end
else;QsWslope=A*0;QsWon=A*0;Hs=A*0;wavePERIOD=A*0;Uwave=A*0;waveANGLE=A*0;PWswell=0*A;deltaPW=0*A;end



%SeaWaves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSeaWaves==1
MASK=0*A+1;
MASK(ho<=hwSea_lim | A==0 | VEG==1)=0;
[Uwave_sea,Tp_sea,Hsea,Fetch,kwave,PWsea]=SeaWaves(h,angleWIND,hwSea_lim,Trange,wind,MASK,64,dx); 
Uwave_sea=Uwave_sea.*(VEG==0 & S==0); Hsea=Hsea.*(VEG==0 & S==0); %vegetation effect and no waves in isolated pond 9because we also redcued 
    if computesand==1; QsWslope_sea=WaveSedimentTransport(Hsea,h,kwave,rhos,N,M,Tp_sea,dx,ss,ws1,hwSea_lim,fTide);QsWslope_sea(Hsea==0)=0;end% Only used for downslope
else;Uwave_sea=0*A;Tp_sea=0*A;Hsea=0*A;Fetch=0*A;QsWslope_sea=0*A;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%Wave-induced edge erosion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computeEdgeErosionSea==1 | computeEdgeErosionSwell==1)%
PW=A*0;
    if computeEdgeErosionSea==1
    PW=PW+PWsea*fMFsea.*fTide; %Wave power reduction for hydroperiod
    end
    if computeEdgeErosionSwell==1;
    PW=PW+PWswell.*fMFswell.*fTide;
    end       
[deltaY1,deltaY2,deltaY3,Pedge,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionSTRAT_3sedimentsXXX(PW,z,aw,maxedgeheight,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;%erode the mardsh edge fully

% %Redistribute the eroded sediment 
if computemud==1;EDGESED=diffuseedgesediments((A==1 & Active==1),EdgeERY2,1*ho,dx);Y2=Y2+EDGESED;end
if computesand==1;EDGESED=diffuseedgesediments((A==1 & Active==1),EdgeERY1,1*ho,dx);Y1=Y1+EDGESED;end
else;Pedge=A*0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%CURRENT-DRIVEN TRANSPORT (Tide and River)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computetide==1 | computeriver==1)
    
    U(A==10)=0;Uwave(A==10)=0;Uwave_sea(A==10)=0; %in the river mouth only resuspension from river flow
    
    %(1)Total sediment resupension SAND
    if computesand==1;
            if ebbfloodcorrection==1
                 [E1E]=totalsedimenterosionSANDsine(h,hlimC,Vebb,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
                 [E1F]=totalsedimenterosionSANDsine(h,hlimC,Vflood,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
                 E1=(E1E+E1F)/2;
            else
                 [E1]=totalsedimenterosionSANDsine(h,hlimC,U,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
            end            
            E1(A==2)=0; %needed for b.c.
    end;
    
    %(2)Total sediment resupension MUD
    if computemud==1
            if ebbfloodcorrection==1
                [E2E,E2tideE,E2swellE]=totalsedimenterosionMUDsineBOUNDARY(A,Vebb,MANN,VEG,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,TrangeVEG,computeSeaWaves,computeSwellwave,computeriver);
                [E2F,E2tideF,E2swellF]=totalsedimenterosionMUDsineBOUNDARY(A,Vflood,MANN,VEG,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,TrangeVEG,computeSeaWaves,computeSwellwave,computeriver);
                E2=(E2E+E2F)/2;
                E2tide=(E2tideE+E2tideF)/2;
                E2swell=(E2swellE+E2swellF)/2;
            else
                [E2,E2tide,E2swell]=totalsedimenterosionMUDsineBOUNDARY(A,U,MANN,VEG,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,TrangeVEG,computeSeaWaves,computeSwellwave,computeriver);
            end      
            E2(A==2)=0; %needed for b.c.
          
    %(3)Total sediment resupension ORGANIC
    E3=E2;
    end

    
    %Erosion limiters when more than one sediemtn class is present
    if computesand==1
    E1=E1.*fracY1;%Reduced for fraction of sediemnt
    end
    if computemud==1;    
    E2=E2.*fracY2;    
    E3=E3.*fracY3;   
    end
         
    
    %Advection-Diffusion Sediment transport
    if computesand==1;
        WS=A*0+ws1;
         [EmD1,SSM,FLX1]=sedtran(0,h,A,SPCLcell,0,DiffSsand,h,ho,E1,WS,dx,dt,rbulk1,co1,Ux,Uy,FLX1,fTide,Ttide,URx,URy,UresX,UresY,periodic,computeriver,computetide,residualcurrents,kro,0,0,QsmouthSAND);  
    else;EmD1=0*A;SSM=0*A;end
    SSC1=SSM./h;
    if computemud==1;
        WS=A*0+ws2;
        WS(VEG==1)=wsB;           
        WS(S==1)=ws2;%Just to be sure, it should not be neceesary 
        [EmD2,SSM,FLX2]=sedtran(1,h,A,SPCLcell,DoMUD+10*Hs,DiffSmud,h,ho,E2,WS,dx,dt,rbulk2,co2,Ux,Uy,FLX2,fTide,Ttide,URx,URy,UresX,UresY,periodic,computeriver,computetide,residualcurrents,kro,1,co2mouth*fMFriver,QsmouthMUD);   
    else;EmD2=0*A;SSM=0*A;end
    SSC2=SSM./h;%./fTide;  %the fTide is needed if want to plot C_H instead of C, only needed for intertidal
    if VEGstratigraphy==1;
        WS=A*0+ws2;
        WS(VEG==1)=wsB;
        [EmD3,SSM,FLX3]=sedtran(1,h,A,SPCLcell,DoMUD,DiffSmud,h,ho,E3,WS,dx,dt,rbulk2,co3,Ux,Uy,FLX3,fTide,Ttide,URx,URy,periodic,computeriver,computetide,kro,1,0,0);
    else;EmD3=0*A;SSM=0*A;end
    SSC3=SSM./h;%./fTide;  %the fTide is needed if want to plot C_H instead of C, only needed for intertidal
    
else;EmD1=0*A;Qs1=0*A;tideE1=A*0;E1=A*0;EmD2=0*A;SSC=0*A;SSC1=A*0;SSC2=A*0;SSC3=A*0;end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Bed evolution erosion/depositon from tidal and river transport
z=imposeboundarydepth(A,z,optionBC,NaN);
if computesand==1;    Y1=Y1-dt*EmD1;end
if computemud==1;     Y2=Y2-dt*EmD2;end
if VEGstratigraphy==1;Y3=Y3-dt*EmD3;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Organic accretion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=B.*(S==0);
if AccreteOrganic==1
if VEGstratigraphy==1;
Y3=Y3+B*Korg*dt; %accrete the organic
else; %put it with mud or sand
    if VEGonsand==0
    Y2=Y2+B*Korg*dt; % putorganic on mud!!!
    else
    Y1=Y1+B*Korg*dt; % putorganic on sand!!!
    end
end
KBTOT=KBTOT+sum(B(A==1))*Korg*dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%BED EVOLUTION VERTICAL FLUXES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the  bed using: 1)Current transport 2)Edge Erosion 3)Organic growth
%This will allow to next compute the evolution by divergence!!!
z=zs+(Y1+Y2+Y3);
znew=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%BED EVOLUTION DIVERGENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EVOLUTION OF Y1
if computesand==1;
    Qs1=E1./(ws1*3600*24).*(U+UR*0.3).*max(hlimC,h); %kg/m/s
    [znew,FQsW_L,FQsW_R,longshore,angleshore]=bedevolutionDIVlongshore(fMFswell*deltaPW,U,fTide,A,AW,z,Zlev,wlo,ho,fracY1,N,M,dt,dx,Trange,Qs1,rbulk1,alphaSAND,NaN,computeSwellwave,(fMFswell*QsWslope+fMFsea*QsWslope_sea*downslopeSANDseawaves),fMFswell*QsWon,angleSWELL,waveANGLE,Active,periodic,optionBC,gridDIR,FQsW_L,FQsW_R);
        deltaY1=z-znew;
        Y1=Y1-deltaY1;
else;longshore=0;angleshore=0;end

%EVOLUTION OF Y2 
if computemud==1;
    z=znew;  %NEED TO RE-UDPATE IT
    Qs2=E2tide./(ws2*3600*24).*(U+UR*0.3).*max(0,h); %kg/m/s %[hcorrected]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);  %VEG=(z-msl)>dBlo;
    [znew]=bedcreepponds(z,A,Active,fracY2,crMUD,crMARSH,dx,dt,VEG,S,Qs2,rbulk2,alphaMUD);  %MUD CREEP  MARSH
        deltaY2=z-znew;
        deltaY2(A==2)=0;%DO NOT UPDATE THE BOUNDARY
        Y2=Y2-deltaY2;
end
 
%EVOLUTION OF Y3 
if computemud==1;
    if VEGstratigraphy==1;
    z=znew;  %NEED TO RE-UDPATE IT
    znew=bedcreepponds(z,A,Active,fracY3,crMUD,crMARSH,dx,dt,VEG,S,Qs2,rbulk2,alphaMUD);  %MUD CREEP  MARSH
        deltaY3=z-znew;
        Y3=Y3-deltaY3;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%IMPOSE "NEUMAN" boundary condition for morphodynamics%%%%%%%%%%%%%%%%%%
%Traslate the first boundary cell bed elevetion
if imposeseaboundarydepthmorphoALL==1;
[plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevationALLBOUNDARY(A,plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
else
    %OPTION 1: IF YOU DO NOTHING YOU LEAVE THE BED AT THE SAME ABSOLUTE LEVEL AND DOES NOT CHANGE WITH RSLR
    %OPTION 2: The bed elevation changes with SLR, so you keep the same depth but not the same bed level
    %     if computesand==1
    %         Y1(A==2)=Y1(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    %     else
    %         Y2(A==2)=Y2(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    %     end
end


%%%%%%%%%%%%%%%%%END OF MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%% NUMERICAL CHECKS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaz=znew-zoriginal;
maxdeltaz=max(abs(deltaz(:)));

muchup=max(0,max(0,znew-zoriginal)).*((znew)>(msl+Trange/2+Hsurge));%modieif on Oct 2019
maxup=max(muchup(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%make an additional smoothing when changes are too large
DIFFsmooth=2;
a=find(abs(deltaz)>limitdeltaz*1 | muchup>limitmaxup*1);
s1=length(a);
if length(a)>1;
        Ytot=(max(0,Y1)+max(0,Y2));fracY1=max(0,Y1)./Ytot;fracY1(Ytot==0)=0;
        if computemud==1
            fracY1(Y1<=0 & Y2<=0)=0.5;
        else
            fracY1=A*0+1;
        end
    muchupY1=A*0;muchupY1(a)=(deltaz(a)).*fracY1(a);
    muchupY2=A*0;muchupY2(a)=(deltaz(a)).*(1-fracY1(a));
Y1=Y1-muchupY1;  DISTRmuchup=diffuseedgesediments((A==1 & Active==1),muchupY1,DIFFsmooth,dx); Y1=Y1+DISTRmuchup;
if computemud==1
Y2=Y2-muchupY2;  DISTRmuchup=diffuseedgesediments((A==1 & Active==1),muchupY2,DIFFsmooth,dx); Y2=Y2+DISTRmuchup;
end
end


znew=zs+(Y1+Y2+Y3);

deltaz=znew-zoriginal;
maxdeltaz=max(abs(deltaz(:)));

muchup=max(0,max(0,znew-zoriginal)).*((znew)>(msl+Trange/2+Hsurge));%modieif on Oct 2019
maxup=max(muchup(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%
if evolvestratigraphy==1;
[Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr]=stratigraphy2D_3sediments(A,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr,nlyr,dlyr,tlyrU,tlyrD);
end
%%%




names = fieldnames(IO);
for i=1:length(names);eval(['IO.' names{i} '=' names{i} ';' ]);end

names = fieldnames(fIO);
for i=1:length(names);eval(['fIO.' names{i} '=' names{i} ';' ]);end

%PLT.PW=PW;
%PLT.PWbank=PWbank;
PLT.angleshore=angleshore;
PLT.wl=wl;
PLT.wlo=wlo;
%PLT.hwaverunup=hwaverunup;
if calculateponddynamics==1
PLT.DIF=DIF;
PLT.DIFo=DIFo;
PLT.dtide=dtide;
PLT.dtideI=dtideI;
end
%PLT.Fetch=Fetch;
%PLOT outputs
PLT.U=U;
PLT.Vebb=Vebb;
PLT.Vflood=Vflood;
PLT.UR=UR;
%if computetide==1 | computeriver==1;
PLT.SSC1=SSC1;
PLT.SSC2=SSC2;
PLT.SSC3=SSC3;
%end
PLT.Tp_sea=Tp_sea;
PLT.Uwave_sea=Uwave_sea;
PLT.Hsea=Hsea;
PLT.Fetch=Fetch;
PLT.EmD2=EmD2;
PLT.Pedge=Pedge;
%PLT.SALTC=SALTC;
%PLT.QsWon=QsWon;
PLT.URx=URx;
PLT.URy=URy;
PLT.Ux=Ux;
PLT.Uy=Uy;
    PLT.Hs=Hs;
if computeSwellwave==1;
    PLT.waveANGLE=waveANGLE;
    PLT.wavePERIOD=wavePERIOD;
    %PLT.kwave=kwave;
    PLT.Uwave=Uwave;
end
if ebbfloodcorrection==1
%PLT.Vx=Vx;PLT.Vy=Vy;
    PLT.Vebb=Vebb;PLT.Vflood=Vflood;
    PLT.VebbX=VebbX;PLT.VfloodX=VfloodX;
    PLT.VebbY=VebbY;PLT.VfloodY=VfloodY;
end
PLT.VEG=VEG;
%PLT.MARSH=MARSH;
%PLT.Ecurv=Ecurv;
%PLT.Cmin=Cmin;
%PLT.Cmax=Cmax;
PLT.deltaPW=deltaPW;
%PLT.HCURV=HCURV;
%PLT.CD=CD;
PLT.wsB=wsB;
%PLT.hriver=hriver;
PLT.h=h;
PLT.AC=AC;
PLT.S=S;
PLT.B=B;
PLT.fTide=fTide;
%PLT.fTidePOND=fTidePOND;
%[angleSWELL ]
PLT.EmD1=EmD1;
%PLT.hdrain=hdrain;
%PLT.qswell=qswell;
%PLT.E1=E1;
PLT.hpRIV=hpRIV;
PLT.ho=ho;
PLT.longshore=longshore;


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            