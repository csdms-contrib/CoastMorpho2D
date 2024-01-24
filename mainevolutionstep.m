function [IO,fIO,maxdeltaz,maxup,PLT,MMM]=mainevolutionstep(A,AW,SPCLcell,PARAMS,x,dx,dt,IO,fIO)

%Load the parameters and variables
names = fieldnames(PARAMS);for i=1:length(names);eval([names{i} '=PARAMS.' names{i} ';' ]);end
names = fieldnames(IO);for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end
names = fieldnames(fIO);for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(A);
%optionBC=2; %impose water depth (for sea boundary)

Y=Y1+Y2+Y3;%total thickness of active later
zs=-zb+(Yb+plyr.*(dlyr)); %the heigth of the ground just below Y
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
    






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE-Calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sea level%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msl=msl+RSLR*dt;
Zlev=z-msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Accrete in the river mouth to maintain same water depth%%%%%%%%%%%%%%%%%%%
if AcrreteRivermouthwithRSLR==1
Y1(A>=10 & A<=19)=Y1(A>=10 & A<=19)+RSLR*dt; %sediment on the bed of the river mouth!!
z=zs+(Y1+Y2+Y3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tidal range attenuation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tidalrangeattenuation==1
[Trange]=gettidalrange(Ttide,Trange_o,CbT1D,Cv,dx,x,Zlev,tempdeltaMSL,ANGLEtideprop,extrapadd);
else 
Trange=A*0+Trange_o;
end
PLT.Trange=Trange;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vegeation upper%limit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dynamicvegationrange==1
    if tidalrangeattenuation==1
    [Trange90]=gettidalrange(Ttide90,Trange90_o,CbT1D,Cv,dx,x,Zlev,tempdeltaMSL,ANGLEtideprop,extrapadd);
    else
    Trange90=Trange90_o;
    end
dBup=MSL90+Trange90/2+0.1;
PLT.Trange90=Trange90;
PLT.dBup=dBup;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pond dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need to run it first because it affects the vegetation through S !!!!!
if calculateponddynamics==1     %%CHECK THE SIGN OF z!!!
zsill=NaN;%TrangeVEG/2*NaN;  
%use lev instead of z to be rleative to MSL!!!
Apond=A;if variablePONDdynamic==1;Apond(pondVAR==0)=0;end;
[deltaY2,pondloss]=pondformation(Apond,dx,dt,Epondform,z-msl,zpondcr,maxdpond,zsill,pondloss,Active,B);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
minponddepthVAR=min(minponddepth,Trange/2);
[S,AC,DIF]=findisolatedponds(z-msl,Apond,N,M,dx,zntwrk,zsill,distdr,minponddepthVAR,Active);%AC is only used for plotting, not in the computation  %DIF is the amount of water impunded at low water. It is the remainng water depth in the pond!
%%%BAA$DS(bio<=0)=0;%you cannot have ponds if it is not vegeated!!!
%S(lev<dBlo)=0;% the ponds in the mudflats are not really a pond! You need to put it otherwise probelm with bedcreeppond. YOU CANNOT CREEP IN PONDS
[S,deltaY2,pondloss]=isolatedpondexpansion(z-msl,S,Apond,N,M,dx,dt,zpondcr,maxdpond,aPEXP,pondloss,Active,B);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
[deltaY2,pondloss]=isolatedponddeepening(z-msl,S,ponddeeprate,dt,pondloss,dBlo);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
else;AC=A*0;S=A*0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vegetation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if VEGETATION==1;
Beq=4*(Zlev-dBup).*(dBlo-Zlev)./(dBlo-dBup).^2;Beq(Zlev<dBlo)=0;Beq(Zlev>dBup)=0;%0.001;  %B=(dBup-lev)/(dBup-dBlo);B(lev>dBup)=0;B(lev<dBlo)=0;
Beq=Beq.*(S==0);%no biomass where there are ponds
%VEG=VEG.*(S==0);%no vegeation where there are ponds %%%OCIOO!!! In some cases I did NOT put theVEG=VEG.*(S==0), so that the
%%%ponds count as vegeated area below (i.e., the sediment transport, and %%%maybe the creep). Don't remember why. Maybe for stability? Check lower%%%tidal range???
PLT.Beq=Beq;

if Dynamicvegeationgrowth==1
[B]=plantbiomass_seeding(A,dx,dt,bioseed,z,B,Zlev,dBseed);
[B]=plantbiomass_lateralexpansion(z,S,A,N,M,dx,dt,plantexpansionrate,B,Active);
B=B./(1-(dt*abiogrow*(1-B./Beq)));% % +dt*abiogrow*Beq)./(1+dt*abiogrow);
B(isnan(B))=0;
else% Static vegeation - only fucntion of elevation
B=Beq;
end

%presence-absence of vegeation
VEG=((B>0 | Zlev>dBup) & A==1 & Active==1);
else;VEG=A*0;Beq=A*0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Keep track of marsh age - does not affect%computation%%%%%%%%%%%%%%%%%%%%%
if trackmarshage==1;
AGE=AGE+dt/365;
AGE(Zlev<dBlo)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Manning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MannH=0*A+CbT;
MannH(VEG==1)=Cv;

MannHR=0*A+CbR;
MannHR(VEG==1)=Cv;

MannSmud=0*A+CbS_MUD;
MannSmud(VEG==1)=CbSveg_MUD;
MannSsand=0*A+CbS_SAND;
MannSsand(VEG==1)=CbSveg_SAND;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%HYDRODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Water depth%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeRiver==1 & riverwaterlevel==1
%hpRIV=max(0,Zlev+initialhforriverflow);%+5; %at the beginning, just a small water level from the river to avoid extra slopes
if usepreviousFLOW==1 & length(hpRIVo)>1;hpRIV=hpRIVo;Uref=URo;
else
hpRIV=max(0,Zlev+initialhforriverflow);%+5; %at the beginning, just a small water level from the river to avoid extra slopes
end
else;hpRIV=0;end
[h,ho,fTide,DHE,wl,wlo]=getwaterdepth(Trange,msl,z,kro,hpRIV,tempdeltaMSL);
%Water depth%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%River flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if computeRiver==1;
    URo=A*0+Uref;%first attempt of velocity
    %URo=FcrUR*sqrt(9.81*h)*0.5;%first attempt of velocity
    DD1=A*0;DDN=A*0;
    [UR,URx,URy,hriver,qR1,qRm1,qRN,qRmN,UR1,URm1,URN,URmN]=riverFLOWiter(A,MannHR,h,ho,dx,fTide,Qmouth,URo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit);%hriver=max(hriver, max(0,-z));

        if riverwaterlevel==1 | correctmomentumriver==1
            URo=UR;
            %hpRIVold=hpRIV;
            %hpRIV=(hpRIV+hriver)/2;
            if riverwaterlevel==0;nitero=1;end
            
                       
%             for i=1:nitero               
%             if riverwaterlevel==1  
%             hpRIV=(hpRIV+hriver)/2;
%             [h,ho,fTide,DHE,wl,wlo]=getwaterdepth(Trange,msl,z,kro,hpRIV,tempdeltaMSL);  
%             if nonlinearfriction==1;URo=max(0.1,(URo*0.9+UR*0.1));end
%             end
%             if  correctmomentumriver==1 & i==nitero; %only to it the last iteration
%             DD1=facmomentum*(UR1-URm1)/dx.*(h>0); DDN=facmomentum*(URN-URmN)/dx.*(h>0);
%             end
%             
%             [UR,URx,URy,hriver,qR1,qRm1,qRN,qRmN,UR1,URm1,URN,URmN]=riverFLOWiter(A,MannHR,h,ho,dx,fTide,Qmouth,URo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit);%hriver=max(hriver, max(0,-z));PLT.hriver=hriver;
%             end            
%             %EX=max((h(:)>1).*abs(hriver(:)-hpRIV(:)))
        
            EX=999;cont=0;%href=h;
            while EX>errorWL & cont<50           
            if riverwaterlevel==1  
            fac=max(factorRiverIter,min(h/10,0.5));%          
            hpRIV=fac.*hriver  +(1-fac).*hpRIV;%  min(0.5,abs(hriver-hpRIV)/2).*sign(hriver-hpRIV);%(hpRIV*0.8+hriver*0.2);
            %hpRIV=hpRIV+ (hriver-hpRIV)/2*min(1,1./max(abs(hriver(:)-hpRIV(:))));%(hpRIV*0.8+hriver*0.2);
            [h,ho,fTide,DHE,wl,wlo]=getwaterdepth(Trange,msl,z,kro,hpRIV,tempdeltaMSL);  
            %if cont==0;href=h;end
            
            %minUfricACT=max(minUfric,0.5*(h<1));
            %if nonlinearfriction==1;URo=max(minUfricACT,(URo*0.5+UR*0.5));end
            if nonlinearfriction==1;URo=max(minUfric,(URo*0.5+UR*0.5));end
            end
            %if  correctmomentumriver==1 & i==nitero; %only to it the last iteration
            %DD1=facmomentum*(UR1-URm1)/dx.*(h>0); DDN=facmomentum*(URN-URmN)/dx.*(h>0);
            %end
            
            [UR,URx,URy,hriver,qR1,qRm1,qRN,qRmN,UR1,URm1,URN,URmN]=riverFLOWiter(A,MannHR,h,ho,dx,fTide,Qmouth,URo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit);%hriver=max(hriver, max(0,-z));PLT.hriver=hriver;
            
            EX=max((h(:)>checkerrorh).*abs(hriver(:)-hpRIV(:)));            %EX= prctile((h(:)>2 & A(:)~=0).*abs(hriver(:)-hpRIV(:)),99.99);
            EX=max(EX,errorUUfac*max((UR(:)>checkerrorU).*(h(:)>checkerrorh).*abs(URo(:)-UR(:))));            %EX= prctile((h(:)>2 & A(:)~=0).*abs(hriver(:)-hpRIV(:)),99.99);            
            cont=cont+1;  
            end         
            [cont EX max(UR(:))]%figure;imagesc(hriver-hpRIV);pause
            
            
%             if  correctmomentumriver==1
%             hlimitDDU=1;
%             DD1=facmomentumRIVER*(UR1-URm1)/dx.*(h>hlimitDDU & A==1); DDN=facmomentumRIVER*(URN-URmN)/dx.*(h>hlimitDDU & A==1);%%%%DD1(A>=10 & A<=19)=0;DDN(A>=10 & A<=19)=0;
%             PLT.DD1=DD1;PLT.DDN=DDN;
%             %if nonlinearfriction==1;URo=max(minUfric,(URo*0.5+UR*0.5));end
%             [UR,URx,URy,hriver,qR1,qRm1,qRN,qRmN,UR1,URm1,URN,URmN]=riverFLOWiter(A,MannHR,h,ho,dx,fTide,Qmouth,URo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit);%hriver=max(hriver, max(0,-z));PLT.hriver=hriver;
%             end

        
        else;hpRIV=A*0;end
        PLT.hriver=hriver;%(hpRIV+hriver)/2;%hriver;
        if usepreviousFLOW==1
        hpRIVo=hpRIV;IO.hpRIVo=hpRIVo;
        IO.URo=URo;
        end
                
        VEGfp=(Zlev>2 & (h<1 | UR<0.5));
        
else;UR=0*A;URy=0*A;URx=0*A;qR1=NaN;qRm1=NaN;qRN=NaN;qRmN=NaN; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if riverwaterlevel==0
% h(A==10)=msl-z(A==10);% the river is not influenced by tide or Hsurge (these two!)
% end

% figure
% subplot(1,2,1);imagesc(URx);colormap('jet');caxis([0 2])
% subplot(1,2,2);imagesc(URx-URxo);colormap('jet');caxis([-0.5 0.5])
% pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-calculate sediment input at river mouth%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeRiver==1;
if computesand
    if imposeSANDCeqartivermouth==1
        for i=1:length(Qmouth)
        [E,Ceq]=totalsedimenterosionSANDsine(hmouth(i),hlimC,kro,CbS_SAND,rho,rhos,ss,d50_1,ws1,1,0,[],[],0,[],[],[],computeRiver,Qmouth(i)/hmouth(i),fMFriver);
        co1mouth(i)=mean(Ceq);%[QsmouthSAND]=Ceq.*h(A==10).*URx(A==10);
        end
    end
else 
end
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the upland cells active
Active(h>kro)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tide&Surge flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetidalcurrent==1;

    if calculateponddynamics==1; %DIF is the impounded water depth, need to subtract to the imput discharge from the ponded area (THE WATER REMAINS THERE!!!)
    DIF=max(DIF,0); %you cannot impound a negative water depth!!! This happens because of the trick used to swap the cell during the pond floodin
    pondleakage=0.2;
    DHE(S==1)=max(0, Trange(S==1)/2-(z(S==1)-msl)-max(0,DIF(S==1)-pondleakage));
    %reduce the hydperperio din the connecte dpond, becuase they do not exahcneg water as much as their wwater depth. Some of that depth is as if it was made of soil. only the top layer count as moving water!
    [~,~,fTidePOND]=getwaterdepth(Trange,msl,z+DIF,kro,hpRIV,tempdeltaMSL);fTide(S==1)=fTidePOND(S==1);%fTide(S==1)=0.01;%
    end

%[U,Uy,Ux,U1,Um1,UN,UmN]=tidalFLOW(A,MannH,h,ho,dx,DHE,Ttide,periodic,kro,0,NaN);
[U,Uy,Ux,U1,Um1,UN,UmN]=tidalFLOWasymmetry(A,MannH,h,ho,dx,DHE,Ttide,periodic,kro,A*0,A*0,DDUlimit,0,NaN);
U(A==10)=0;Ux(A==10)=0;Uy(A==10)=0;
Ubase=U;PLT.Ubase=Ubase;
else;U=0*A;Uy=0*A;Ux=0*A;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shallow tidal%flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calculateshallowflow==1;
[hS,hoS]=getwaterdepthSHALLOW(Trange,msl,z,kroS,hpRIV,tempdeltaMSL);
DHES=(Zlev<dBlo & Zlev>(-Trange/2) & S==0).*Trange/nSHALLOW*pi/2;
%[US]=tidalFLOW(A,MannH,hS,hoS,dx,DHES,Ttide/nSHALLOW,periodic,kroS,0,NaN);
[US]=tidalFLOWasymmetry(A,MannH,hS,hoS,dx,DHES,Ttide/nSHALLOW,periodic,kro,A*0,A*0,DDUlimit,0,NaN);
US=min(maxQshallow,hS.*US)./hS;
US=min(US,FcrUS*sqrt(9.81*hS));
US(VEG==1)=0;
PLT.US=US;PLT.hS=hS;
else
US=A*0;hS=A*0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CORRECTION FOR% CURVATUE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetidalcurrent==1 & curvaturecorrection==1;

MASKCURVATURE=(U>0.1 & VEG==0) ;%;subplot(2,1,2);%[FX,FY] = gradient(z);slope=sqrt(FX.^2+FY.^2)/dx;%MASKCURVATURE=(slope<0.05 & VEG==0);

%Normal velocity vector
Unx=-Uy;Uny=Ux;
Undiffx=-Uy;Undiffy=Ux;

Umagnitude=sqrt(Ux.^2+Uy.^2);
Uxo=Ux./Umagnitude;Uyo=Uy./Umagnitude;
Uxo(MASKCURVATURE==0)=NaN;Uyo(MASKCURVATURE==0)=NaN;

cur=curlNAN(Uyo,Uxo)/dx.*U.*h;
%cur=sqrt(abs(cur)).*sign(cur);
cur=abs(cur).^0.2.*sign(cur);
cur=diffusecurvature(isfinite(cur),cur,a_diffusecur,dx,Unx,Uny);
cur(isnan(cur))=0;

Uno=sqrt(Unx.^2+Uny.^2);Unx=Unx./Uno.*cur;Uny=Uny./Uno.*cur;
Unx(isnan(Unx))=0;Uny(isnan(Uny))=0;

%Factor for flow MODIFICATION -implemented after
factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx,Uny,MASKCURVATURE,Undiffx,Undiffy); 
factorU(MASKCURVATURE==1)=min(1,factorU(MASKCURVATURE==1).^0.5);
factorU(MASKCURVATURE==0)=1;

%STOCASTIC BANK EROSION
if flowbankerosion==1
MASK=MASKCURVATURE;
Pbank=abs(cur)*dx;
a_bankerosion=A*0+a_bankerosion;
a_bankerosion(VEG==0)=a_bankerosionUNVEG;
[deltaY1,deltaY2,deltaY3,PedgeBANK,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionBANKPUSH(Pbank,Unx,Uny,z,a_bankerosion,999,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;%erode the mardsh edge fully
EDGESED=diffuseedgesediments((A==1),EdgeERY2,0.1*ho.*Ubase,dx);Y2=Y2+EDGESED; %put eroded sediment only in chnannels, not on marsh
EDGESED=diffuseedgesediments((A==1),EdgeERY1,0.1*ho.*Ubase,dx);Y1=Y1+EDGESED; %put eroded sediment only in chnannels, not on marsh
end

PLT.cur=cur;PLT.factorU=factorU;  PLT.MASKCURVATURE=MASKCURVATURE;
else
factorU=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CORRECTION FOR 1)Curvature, 2)Ebb-Flood, 3) Nonliner fricton%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetidalcurrent==1 & (correctmomentumtideebbflood==1 | tidalnonlinearflow==1 | curvaturecorrection==1)
    
if correctmomentumtideebbflood==0 %This also consider the non-linear friction if specified
[U,Uy,Ux]=tidalFLOWasymmetry(A,MannH./factorU,h,ho,dx,DHE,Ttide,periodic,kro,A*0,A*0,DDUlimit,tidalnonlinearflow,Ubase);
elseif correctmomentumtideebbflood==1
hlimitDDU=1;
DD1=facmomentumTIDE*(U1-Um1)/dx.*(h>hlimitDDU); DDN=facmomentumTIDE*(UN-UmN)/dx.*(h>hlimitDDU);
DD1(A==2)=0;DDN(A==2)=0;
[Vebb,Vebby,Vebbx]=tidalFLOWasymmetry(A,MannH./factorU,h,ho,dx,DHE,Ttide,periodic,kro,DD1,DDN,DDUlimit,tidalnonlinearflow,Ubase);
[Vflood,Vfloody,Vfloodx]=tidalFLOWasymmetry(A,MannH./factorU,h,ho,dx,DHE,Ttide,periodic,kro,-DD1,-DDN,DDUlimit,tidalnonlinearflow,Ubase);
Ux=(Vebbx+Vfloodx)/2;Uy=(Vebby+Vfloody)/2;
U=sqrt(Ux.^2+Uy.^2);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Connect the channel corners%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if connectchannelcorners==1
MASKflow=double(VEG==0);
MASKflow1=bwmorph(MASKflow,'bridge');
MASKflow2=bwmorph(MASKflow1,'diag');
fU=0.2;
Ufmax=max(U/fU,[U(:,1) U(:,1:end-1)]);
Ufmax=max(Ufmax,[U(:,2:end) U(:,end)]);
Ufmax=max(Ufmax,[U(1,:); U(1:end-1,:)]);
Ufmax=max(Ufmax,[U(2:end,:); U(end,:)]);
U(MASKflow2==1 & MASKflow==0)=fU*Ufmax(MASKflow2==1 & MASKflow==0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CORRECTION FOR%MOMENTUM%ADVECTION$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calculatesalinity==1;
[~,SALTM]=sedtran(1,h,A,SPCLcell,0,DiffSmud,Zlev,h,ho,A*0,A*0,dx,dt,rbulk1,SALTocean,0,NaN,Ux,Uy,fTide,Ttide,qR1,qRm1,qRN,qRmN,periodic,computeRiver,computetidalcurrent,residualcurrents,kro,co1mouth*0,FLX1,0,[]);           
SALTC=SALTM./h;
PLT.SALTC=SALTC;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Swell waves%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSwellWave==1;   
       
    hwave=ho;%just to isolate this water depth and not mess up    
%     if Hsurge>0
%     hwaverunup=(max(1,h)-h);
%     hwave=hwave+hwaverunup;
%     wlo=wlo+hwaverunup;
%     end   
   
        %THIS IS THE OLD SINGLE FREQUENCY AND SINGLE DIRECTION
         % kwave=wavek(1/Tp_swell,hwave);
         % [Hs]=SwellWaves(A,AW,Ho,N,M,hwave,Tp_swell,kwave,dx,periodic,angleSWELL,gridDIR);  
    if multifrequency==1; %Case multiple frequency
        [Tperiodi Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,[1 1.5 2 2.5]);
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbrk,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap,nrefrac,wavediffraction);
    elseif multifrequency==0 %Case single frequency  
         %%%%%%%%%%Tperiodi=Tp_swell;Ejonswap=1;[Hs,waveANGLE,wavePERIOD,PWswell]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,hwave,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap);
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirection(A,AW,Ho,N,M,Cbrk,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tp_swell,dx,periodic,angleSWELL,gridDIR,nrefrac,wavediffraction);
    end
   waveANGLE(isnan(waveANGLE))=0;
   %Hs=0*h;Hs(h>0.5)=Ho;%to reproduce ortiz results
   %figure;
   %subplot(1,2,1);imagesc(PWswell);set(gca,'YDir','normal');colormap('jet');
   %subplot(1,2,2);imagesc(deltaPW);set(gca,'YDir','normal');colormap('jet');caxis([0 20]);
   %pause
          if computesand==1;
           [QsWslope,QsWon]=WaveSedimentTransport(Hs,hwave,kwave,rhos,N,M,wavePERIOD,dx,ss,ws1,hwSwell_lim,fTide);
%            QsWon(hwave<=hwSwelltransport_lim)=0;
%            QsWslope(hwave<=hwSwelltransport_lim)=0;
%            deltaPW(hwave<=hwSwelltransport_lim)=0;   
           QsWon(hwave<=hwSwell_lim)=0;
           QsWslope(hwave<=hwSwell_lim)=0;
           deltaPW(hwave<=hwSwell_lim)=0;
           
           %PWswell(hwave<=hwSwell_lim)=0;
          end
else;QsWslope=A*0;QsWon=A*0;Hs=A*0;wavePERIOD=A*0;Uwave=A*0;waveANGLE=A*0;PWswell=0*A;deltaPW=A*0;
end
%if periodic==0
%mA=mean(waveANGLE(:));waveANGLE(waveANGLE*mA<1)=0.001*sign(mA);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%swell wave effect on destorying marshes
if VEGETATION==1 & computeSwellWave & swelldestroyVEG==1
B(Hs>0.5)=0;
VEG(Hs>0.5)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSeaWave==1
%hwave=h;   %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
%MASK=0*A+1;MASK(hwave<=hwSea_lim | A==0 | VEG==1)=0;
%[Uwave_seaI,Tp_seaI,Hsea,Fetch,PWsea,QsWslope_sea]=SeaWaves_multipleheight(hwave,WINDdir,hwSea_lim,Trange,WINDspeed,MASK,64,dx,z,msl,tempdeltaMSL,ws1,fTide,extraHseaimposed,addextrafetch,extrafetch); %72,dx
[Uwave_seaI,Tp_seaI,Hsea,Fetch,PWsea,QsWslope_sea]=SeaWaves_multipleheight(A,WINDdir,hwSea_lim,Trange,WINDspeed,VEG,64,Nhseawave,dx,z,msl,tempdeltaMSL,ws1,fTide,extraHseaimposed,addextrafetch,extrafetch); %72,dx

Uwave_seaI=Uwave_seaI.*(VEG==0 & S==0); Hsea=Hsea.*(VEG==0 & S==0); %vegetation effect and no waves in isolated pond 9because we also redcued ws!!1)%Uwave_sea=Uwave_sea.*(VEG==0); Hsea=Hsea.*(VEG==0); %vegetation effect
else;Uwave_sea=0*A;Hsea=0*A;Fetch=0*A;QsWslope_sea=0*A;Uwave_seaI=NaN;Tp_seaI=NaN;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size(Uwave_seaI)
%figure;imagesc(Uwave_seaI(:,:,1));pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wave-induced edge erosion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computeEdgeErosionSea==1 | computeEdgeErosionSwell==1)%%%MASK=0*A+1;MASK(h<hwSea_lim | A==0)=0;
PW=A*0;
    if computeEdgeErosionSea==1
    PW=PW+PWsea;%*fMFsea.*fTide; %Wave power reduction for hydroperiod
    end
%     if computeEdgeErosionSwell==1;
%     PW=PW+PWswell.*fMFswell.*fTide;
%     end       
%MASK=0*A+1;MASK(A==0 | VEG==1)=0;
MASK=0*A+1;MASK(h<=hwSea_lim | A==0 | VEG==1)=0;

if variableEDGEEROSION==0;awVAR=[];end

[deltaY1,deltaY2,deltaY3,Pedge,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionSTRAT_3sedimentsXXX(PW,z,aw,maxedgeheight,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX,variableEDGEEROSION,awVAR);
Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;%erode the mardsh edge fully
%Redistribute the eroded sediment  
EDGESED=diffuseedgesediments((A==1),EdgeERY2,0.1*ho,dx);Y2=Y2+EDGESED;
EDGESED=diffuseedgesediments((A==1),EdgeERY1,0.1*ho,dx);Y1=Y1+EDGESED;
%EDGESED=diffuseedgesediments((A==1),EdgeERY2,h,dx);Y2=Y2+EDGESED;
%EDGESED=diffuseedgesediments((A==1),EdgeERY1,h+0.2,dx);Y1=Y1+EDGESED;
else;Pedge=A*0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CURRENT-DRIVEN TRANSPORT (Tide and River)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computetidalcurrent==1 | computeRiver==1)
    
    U(A==10)=0;Uwave(A==10)=0;Uwave_sea(A==10)=0; %in the river mouth only resuspension from river flow
    
    %(1)Total sediment resupension SAND
    if computesand==1;
            if correctmomentumtideebbflood==1
                 [E1E]=totalsedimenterosionSANDsine(h,hlimC,kro,MannSsand,rho,rhos,ss,d50_1,ws1,fTide,computetidalcurrent,Vebb,FcrUT,calculateshallowflow,US,hS,nSHALLOW,computeRiver,UR,fMFriver);
                 [E1F]=totalsedimenterosionSANDsine(h,hlimC,kro,MannSsand,rho,rhos,ss,d50_1,ws1,fTide,computetidalcurrent,Vflood,FcrUT,calculateshallowflow,US,hS,nSHALLOW,computeRiver,UR,fMFriver);
                 E1=(E1E+E1F)/2;
            else
                 [E1]=totalsedimenterosionSANDsine(h,hlimC,kro,MannSsand,rho,rhos,ss,d50_1,ws1,fTide,computetidalcurrent,U,FcrUT,calculateshallowflow,US,hS,nSHALLOW,computeRiver,UR,fMFriver);
            end            
            E1(A==2)=0; %needed for b.c.  
            PLT.E1=E1;
            %E1=min(E1/rbulk1*dt,5)*(rbulk1/dt);
    end;
    
    %(2)Total sediment resupension MUD
    if computemud==1
                    if computeSwellWave==1;%Swell waves for MUD only!
                    UwaveMUD=facHwave*Uwave.*(VEG==0); %vegetation effect. Plants put to zero wave erosion
                    else;QsWslope=zeros(N,M);QsWon=zeros(N,M);BRK=zeros(N,M);UwaveMUD=zeros(N,M);Hs=zeros(N,M);Tp_swellMUD=1;HswellMUD=0;end
            if correctmomentumtideebbflood==1
                  [E2E,~,B]=totalsedimenterosionMUDsine(ho,h,taucr,taucrVEG,VEG,me,MannSmud,fTide,Vebb,FcrUT,FcrUTveg,hlimC,calculateshallowflow,US,hS,nSHALLOW,computeSeaWave,Uwave_seaI,Tp_seaI,computeSwellWave,UwaveMUD,Tp_swell,computeRiver,fMFriver,UR,flowdestroyVEG,B);
                  [E2F,~,B]=totalsedimenterosionMUDsine(ho,h,taucr,taucrVEG,VEG,me,MannSmud,fTide,Vflood,FcrUT,FcrUTveg,hlimC,calculateshallowflow,US,hS,nSHALLOW,computeSeaWave,Uwave_seaI,Tp_seaI,computeSwellWave,UwaveMUD,Tp_swell,computeRiver,fMFriver,UR,flowdestroyVEG,B);
                  E2=(E2E+E2F)/2;
            else          
                 [E2,~,B]=totalsedimenterosionMUDsine(ho,h,taucr,taucrVEG,VEG,me,MannSmud,fTide,U,FcrUT,FcrUTveg,hlimC,calculateshallowflow,US,hS,nSHALLOW,computeSeaWave,Uwave_seaI,Tp_seaI,computeSwellWave,UwaveMUD,Tp_swell,computeRiver,fMFriver,UR,flowdestroyVEG,B);
            end      
            E2(A==2)=0; %needed for b.c.      
            PLT.E2=E2;
          
    %(3)Total sediment resupension ORGANIC
    E3=E2;
    end
    

    
%Erosion limiters
    if computesand==1
    E1o=E1;
    fracY1(A>=10 | A<=19)=1;
    E1=E1.*fracY1;%Reduced for fraction of sediemnt due to coesion
    Elimit=max(0,Y1*conces)/dt*rbulk1;a=find(E1>Elimit & A==1);E1(a)=Elimit(a);%this is the limit to avoid E to scour more than Y1 or Y2. 
    end
    if computemud==1;  
    E2o=E2;  
    % conces=1;%how much to extra erode, a parameter
    %E2=E2.*fracY2;    
    Elimit=max(0,Y2*conces)/dt*rbulk2;a=find(E2>Elimit);E2(a)=Elimit(a);
    %E2(Y2<0.2 & VEG==1)=0;
    %E3=E3.*fracY3;   
    %Elimit=max(0,Y3*conces)/dt*rbulk2;a=find(E3>Elimit);E3(a)=Elimit(a);
    end
    

%Advection-Diffusion Sediment transport
    if computesand==1;
        WS=A*0+ws1;
        flagSANDMUD=0;
           [EmD1,SSM,FLX1,XXTide,XXRiver]=sedtran(flagSANDMUD,h,A,SPCLcell,DoSAND,DiffSsand,Zlev,h,ho,E1,WS,dx,dt,rbulk1,co1,SeaSSCbelowLIMIT,ZlevcoLIMIT,Ux,Uy,fTide,Ttide,qR1,qRm1,qRN,qRmN,periodic,computeRiver,computetidalcurrent,residualcurrents,kro,co1mouth,FLX1,tracksedimentfluxes,XX); 
           %[EmD1,SSM,FLX1,XXTide,XXRiver]=sedtran(flagSANDMUD,h,A,SPCLcell,DoMUD,DiffSsand,Zlev,h,ho,E1,WS,dx,dt,rbulk1,co1,SeaSSCbelowLIMIT,ZlevcoLIMIT,Ux,Uy,fTide,Ttide,qR1,qRm1,qRN,qRmN,periodic,computeRiver,computetidalcurrent,residualcurrents,kro,co1mouth,FLX1,tracksedimentfluxes,XX); 
           PLT.XXTide1=XXTide;PLT.XXRiver1=XXRiver;
           if (computeRiver==1 & imposeNOerosiondepostionatmouthSAND==1);EmD1(A>=10 & A<=19)=0;end
    else;EmD1=0*A;SSM=0*A;end
    SSC1=SSM./h;
    if computemud==1;
        WS=A*0+ws2;
        WS(VEG==1)=wsB;           
        WS(S==1)=ws2;%SHOULD NOT BE NECEEARY BECUASE VEG alreeady set equal to zero where S=1 (see above).  ->Do not add the vegetation settling velocity in the ponds! %WS(S==1)=0.000000000001;%E2(S==1)=0;
        flagSANDMUD=1;  
        ADoMUD=A*0+DoMUD;
        ADoMUD(VEG==1)=DoMUDveg;
        ADoMUD(VEG==1)=ADoMUD(VEG==1)+DoMUDsubgridVEG;
        %PLT.ADoMUD=ADoMUD;
        [EmD2,SSM,FLX2,XXTide,XXRiver]=sedtran(flagSANDMUD,h,A,SPCLcell,ADoMUD,DiffSmud,Zlev,h,ho,E2,WS,dx,dt,rbulk2,co2,0,NaN,Ux,Uy,fTide,Ttide,qR1,qRm1,qRN,qRmN,periodic,computeRiver,computetidalcurrent,residualcurrents,kro,co2mouth,FLX2,tracksedimentfluxes,XX); 
        PLT.XXTide2=XXTide;PLT.XXRiver2=XXRiver;
        if (computeRiver==1 & imposeNOerosiondepostionatmouthMUD==1);EmD2(A>=10 & A<=19)=0;end    
    else;EmD2=0*A;SSM=0*A;end
    SSC2=SSM./h;%./fTide;  %devi metter i fTide per farti plottare la b.c quando il fondo e' sopra il MLW (ftide<1)
    if VEGstratigraphy==1;
        %WS=A*0+ws2;WS(VEG==1)=wsB;
        %flagSANDMUD=1;
        %[EmD3,SSM,FLX3]=sedtran(flagSANDMUD,h,A,SPCLcell,DoMUD,DiffSmud,h,ho,E3,WS,dx,dt,rbulk2,co3,Ux,Uy,FLX3,fTide,Ttide,URx,URy,periodic,computeRiver,computetidalcurrent,kro,1,0,0);
    else;EmD3=0*A;SSM=0*A;end
    SSC3=SSM./h;%./fTide;  %devi metter il fTide per farti plottare la b.c quando il fondo e' sopra il MLW (ftide<1)
    
else;EmD1=0*A;Qs1=0*A;tideE1=A*0;E1=A*0;EmD2=0*A;SSC=0*A;SSC1=A*0;SSC2=A*0;SSC3=A*0;end 
%end of compute tide and river%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Bed evolution erosion/depositon from tidal and river transport
%z=imposeboundarydepth(A,z,optionBC,NaN);
if computesand==1;    Y1=Y1-dt*EmD1;end
if computemud==1;     Y2=Y2-dt*EmD2;end
if VEGstratigraphy==1;Y3=Y3-dt*EmD3;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Impose depth of river
% z=zs+(Y1+Y2+Y3);
% hR=hpRIV(A==10)-z(A==10);
% Y1(A==10)=Y1(A==10)+hR-hmouth(1);
% %%%%




%Organic accretion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=B.*(S==0).*Active.*(A==1);
if AccreteOrganic==1
if variableORGaccretion==1;Korg=KorgVAR;end
if VEGstratigraphy==1
Y3=Y3+B.*Korg*dt; %accrete the organic
else; %put it with mud or sand
    if VEGonsand==0
    Y2=Y2+B.*Korg*dt; % putorganic on mud!!!
    else
    Y1=Y1+B.*Korg*dt; % putorganic on sand!!!
    end
end
BK=B.*Korg;
KBTOT=KBTOT+sum(BK(A==1))*dt;
end


% %%%MARSH COMPACTION
% Y2(VEG==1 & A==1)=Y2(VEG==1 & A==1)-COMP*dt;
% KBTOT=KBTOT-sum(VEG(:)==1 & A(:)==1)*COMP*dt;
Y2(B>0 & A==1)=Y2(B>0 & A==1)-COMP*dt;
%KBTOT=KBTOT-sum(B(:)>0 & A(:)==1)*COMP*dt;
Y2OX=Y2OX+sum(B(:)>0 & A(:)==1)*COMP*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%BED EVOLUTION VERTICAL FLUXES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the  bed using: 1)Current transport 2)Edge Erosion 3)Organic growth
z=zs+(Y1+Y2+Y3);
znew=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%Aeolian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Procesaeolian==1;    
Zlev(1:2,:)=10;
    zaeo1=Trange_o/2;%NON CAMBIARLO 1) non migliora la barrier, 2) fa' peggio nel mainland
    zaeo2=Trange_o/2+3;
    
    %Zdw=[Zlev(1,:); Zlev(1:end-1,:)];
    BB=zeros(M,1);
    for i=1:M;a=find(Zlev(:,i)>zaeo1);
        BB(i)=a(end);
        if Hs(a(end)+floor(500/dx))<0.1;BB(i)=NaN;end
        %Hfront=Hs(a(end)+1);
    end
    bdist=(ones(N,1)*BB'-[1:N]'*ones(1,M));
   
%     %option1
%      %QsAEO=0.05*(10+bdist*dx).^-0.3;%distycee IBL
%      QsAEO=0.02*(1+bdist*dx).^-0.3;%distycee IBL
%      %QsAEO=0.1*(10+bdist*dx).^-0.3;%distycee IBL
%      QsAEO(bdist<0)=0;
%      QsAEO=QsAEO.*(Active==1).*(Zlev>zaeo1);%do not move wet sand
%      QsAEO=QsAEO.*(Zlev>zaeo1 & Zlev<zaeo2).*((zaeo2-Zlev)/(zaeo2-zaeo1)); %so that the tall dunes are not eroded by wind
%      QsAEO=QsAEO.*(Zdw<zaeo2);%do not move if the next cell is too high
     
     %option1b
     ZlevCUM=cummax(Zlev,'reverse');
     ZlevCUM(Zlev<zaeo1)=Zlev(Zlev<zaeo1);%so that the cummax does not "raise" the points lower than the low-limit
     QsAEO=1*0.05*(10+bdist*dx).^-0.3;%distycee IBL
     QsAEO(bdist<0 | isnan(bdist))=0;    
     QsAEO=QsAEO.*(Active==1);%do not move wet sand
     a=find(Active==0);QsAEO(a+1)=0;%do not move if inactive the landward cell too
     QsAEO=QsAEO.*(ZlevCUM>zaeo1 & ZlevCUM<zaeo2).*((zaeo2-ZlevCUM)/(zaeo2-zaeo1)); %so that the tall dunes are not eroded by wind
     
%      Zdw=[Zlev(1,:); Zlev(1:end-1,:)];
%      Zdw=max(0,Zdw-z
%      QsAEO=
     
%     %option 2
%      QsAEO=0.01*(Zlev>zaeo1 & Zlev<zaeo2).*((zaeo2-Zlev)/(zaeo2-zaeo1)); %so that the tall dunes are not eroded by wind
%      QsAEO=QsAEO.*(Active==1);%do not move wet sand
%      QsAEO=QsAEO.*(Zdw<zaeo2);%do not move if the next cell is too high
%      

     
    zoRUNP=0-0.5;
    zmaxRUNP=Trange_o/2+0.8;
        Qrunp=A*0;for i=1:M;a=find(Zlev(:,i)>zoRUNP);Qrunp(a(end),i)=1;end
        a=find(Qrunp==1);
        Qrunp(a(Zlev(a)>zmaxRUNP))=0;%too high to be eroded
        Hrunlim=0.2;Qrunp(a(Hs(a+1)<Hrunlim & Hs(a+2)<Hrunlim & Hs(a+3)<Hrunlim))=0;%not enough wave in front
        Qrunp(a(Zlev(a-1)>zmaxRUNP))=0;%too step to move upward
        %Qrunp=Qrunp*0.1;
        Qrunp=Qrunp*0.1*2;
        PLT.Qrunp=Qrunp;

  %QsDRY=QsAEO+Qrunp;
  QsDRY=Qrunp;
%%%

    deltaA=(QsDRY-[QsDRY(2:end,:); QsDRY(end,:)])/dx*dt;
    Y1=Y1-deltaA;
    znew=znew-deltaA;
    z=znew;
     
%     %diffusion of aeolian
%     dh=max(0,z-msl-zaeo1);
%     dhNEW=diffuseedgesediments((A==1),dh,0.0000005*1*10*(Active==1 & dh>0)*dt,dx);
%     delta=dh-dhNEW;
%     Y1=Y1-delta;
%     znew=znew-delta;
%     z=znew;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%BED EVOLUTION DIVERGENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EVOLUTION OF Y1
if computesand==1;
Yreduction=(1-min(1,exp(-10*(Y1)))).*fracY1;    %ONLY FOR THE WAVE TRANSPORT, DIOCANE
Qs1=E1o./(ws1*3600*24).*(U*alphaSAND+UR*alphaSANDriver).*max(hlimCdwn,h); %kg/m/s% USED FOR BARNSTABLE

%Qs1=E1o./(ws1*3600*24).*(U*alphaSAND+UR*alphaSANDriver).*h.*(h>0.5); %kg/m/s
%Qs1=E1o./(ws1*3600*24).*((Vebb+Vflood)/2*alphaSAND+UR*alphaSANDriver).*max(hlimCdwn,h); %kg/m/s
%wlo is the water level in which the dpeht h can actually go to zer (h is the water level in which the depth is at minimum kro, do avoid explidign with flow and disperive transport. wlo is used in wave transpor wlo=wlo+hwaverunup;%add extra water level due to wave run up. Already added where you calculate the swell waves around line 189
if computeSwellWave==1 | computeSeaWave==1;computewave=1;else;computewave=0;end;
[znew,FQsW_L,FQsW_R,longshore,angleshore]=bedevolutionDIVlongshore(deltaPW,U,fTide,A,AW,z,Zlev,wlo,ho,Y1,Yreduction,VEG,N,M,dt,dx,Trange,Qs1,rbulk1,hwSwelltransport_lim,computewave,QsWslope+QsWslope_sea*downslopeSANDseawaves,QsWon,angleSWELL,waveANGLE,Active,periodic,wavetransportzerolateralgradient,gridDIR,FQsW_L,FQsW_R);
deltaY1=z-znew;Y1=Y1-deltaY1;
else;longshore=0;angleshore=0;end

%EVOLUTION OF Y2 
if computemud==1;
z=znew;  %NEED TO RE-UDPATED IT
Yreduction=(1-min(1,exp(-10*(Y2)))).*fracY2;    %fracY2; %this is just for the marsh stratigraphy%FORSE PRIMA ERA PER UNA NO PER 10    %Yreduction=fracY2;  %this is for most of the simulations big basin % Yreduction(Yreduction<0.5)=0;
Qs2=E2o./(ws2*3600*24).*(U+UR*0.1).*max(hlimCdwnMUD,h); %kg/m/s 

znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S,Qs2,rbulk2,alphaMUD,facQsbank,U+UR);%,deltaUC,a_bankcreep);  %MUD CREEP  MARSH
%znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S,Qs2,rbulk2,alphaMUD,facQsbank,Ubase+UR);%,deltaUC,a_bankcreep);  %MUD CREEP  MARSH
deltaY2=z-znew;
deltaY2(A==2)=0;  %DO NOT UPDATE THE BOUNDARY
Y2=Y2-deltaY2;
end
 
%EVOLUTION OF Y3 
if computemud==1;
    if VEGstratigraphy==1;
    z=znew;  %NEED TO RE-UDPATED FROM THE Y1
    Yreduction=fracY3;%fracY3; (1-min(1,exp(-1*(Y3)))).*%Yreduction(Yreduction<0.1)=0;%to reduce some starneg probelms with creep of two sediment swiht stratigraphy. negative creep!
    znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S);  %MUD CREEP
    deltaY3=z-znew;Y3=Y3-deltaY3;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if evolvestratigraphy==1;
[Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr]=stratigraphy2D_3sediments(A,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr,nlyr,dlyr,tlyrU,tlyrD);
end

%%%%IMPOSE "NEUMAN" boundary condition for morphodynamics%%%%%%%%%%%%%%%%%%
%Traslate the first boundary cell bed elevetion
if imposeseaboundarydepthmorphoALL==1;
[plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,zb]=seaboundaryNeumanbedelevationALLBOUNDARY(A,plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,zb);
else
    %remeber that zb has the oppisite sign!!!
    if computesand==1
        Y1(A==2)=Y1(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    else
        Y2(A==2)=Y2(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    end
end


%%%%%%%%%%%%%%%%%END OF MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%pause


% % %%%%%%%%%%%%%%%%
% znew=z;
% znew(end-50:end,:)=-5+msl;
% deltaY2=z-znew;
% Y2=Y2-deltaY2;
% %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% NUMERICAL
%%%%%%%%%% CHECKS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxdeltaz=prctile(abs(znew(:)-zoriginal(:)),99.9);
%ORGINAL VERSION muchup=max(0,max(0,znew-zoriginal)-DHE);
%MMM=(znew-zoriginal).*((znew)>(msl+Trange/2+tempdeltaMSL));
muchup=max(0,max(0,znew-zoriginal)).*((znew)>(msl+Trange/2+tempdeltaMSL));%modieif on Oct 2019
maxup=max(muchup(:));%maxup=prctile(muchup(:),99.9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
% IO.zb=zb;
% IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
% IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
% IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
% IO.Active=Active;
% fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
% fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
% fIO.pondloss=pondloss;
% fIO.FQsW_R=FQsW_R;
% fIO.FQsW_L=FQsW_L;
names = fieldnames(IO);
for i=1:length(names);eval(['IO.' names{i} '=' names{i} ';' ]);end

names = fieldnames(fIO);
for i=1:length(names);eval(['fIO.' names{i} '=' names{i} ';' ]);end

%PLT.PW=PW;
%PLT.angleshore=angleshore;
PLT.wl=wl;
PLT.wlo=wlo;
%PLT.hwaverunup=hwaverunup;
if calculateponddynamics==1
PLT.DIF=DIF;
%PLT.dtide=dtide;
%PLT.dtideI=dtideI;
end
%PLT.Fetch=Fetch;
%PLOT outputs
PLT.U=U;
PLT.UR=UR;
%if computetidalcurrent==1 | computeRiver==1;
PLT.SSC1=SSC1;
PLT.SSC2=SSC2;
PLT.SSC3=SSC3;
%end
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
if computeSwellWave==1;
    PLT.waveANGLE=waveANGLE;
    PLT.wavePERIOD=wavePERIOD;
    PLT.Uwave=Uwave;
end
if correctmomentumtideebbflood==1
    PLT.Vebb=Vebb;PLT.Vflood=Vflood;
    %PLT.VebbX=VebbX;PLT.VfloodX=VfloodX;
    %PLT.VebbY=VebbY;PLT.VfloodY=VfloodY;
end
PLT.VEG=VEG;
%PLT.MARSH=MARSH;
%PLT.Cmin=Cmin;
%PLT.Cmax=Cmax;
PLT.deltaPW=deltaPW;
%PLT.HCURV=HCURV;
%PLT.CD=CD;
PLT.wsB=wsB;
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
%PLT.longshore=longshore;




% if imposeseaboundarydepthmorphoEAST==1;%east
% [plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevationEAST(plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
% end
% if imposeseaboundarydepthmorphoNORTH==1;%north
% [plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevation(plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
% end
% sumY1=sumSedcolum(Yb,flyrb1,flyr1,dlyr,Y1);sumY1=sum(sumY1(A==1))









% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z=zs-(Y1+Y2+Y3);%update the bed elevation
% [h,ho,fTide,dtide,DHE,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);
% 
% % %Swell waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if computeSwellWave==1;   
%     hwave=ho;%just to isolate this water depth and not mess up
%     
%         %THIS IS THE OLD SINGLE FREQUENCY AND SINGLE DIRECTION
%          % kwave=wavek(1/Tp_swell,hwave);
%          % [Hs]=SwellWaves(A,AW,Ho,N,M,hwave,Tp_swell,kwave,dx,periodic,angleSWELL,gridDIR);   
%     if multifrequency==1; %Case multiple frequency
%         [Tperiodi Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,[1 1.5 2 2.5]);
%         [Hs,waveANGLE,wavePERIOD,PWswell,kwave]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbrk,Cbed,wavefrictionCollins,hwave,hwSwell_lim,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap,nrefrac,wavediffraction);
%     elseif multifrequency==0 %Case single frequency  
%          %%%%%%%%%%Tperiodi=Tp_swell;Ejonswap=1;[Hs,waveANGLE,wavePERIOD,PWswell]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,hwave,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap);
%         [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave]=SwellWavesMultiDirection(A,AW,Ho,N,M,Cbrk,Cbed,wavefrictionCollins,hwave,hwSwell_lim,Tp_swell,dx,periodic,angleSWELL,gridDIR,nrefrac,wavediffraction);
%     end
%    waveANGLE(isnan(waveANGLE))=0;
%    %Hs=0*h;Hs(h>0.5)=Ho;%to reproduce ortiz results
%    
%           if computesand==1;
%            [QsWslope,QsWon]=WaveSedimentTransport(Hs,hwave,kwave,rhos,N,M,wavePERIOD,dx,ss,ws1,hwSwell_lim,fTide);
%            QsWon(hwave<=hwSwell_lim)=0;QsWslope(hwave<=hwSwell_lim)=0;
%           end
% else;QsWslope=A*0;QsWon=A*0;Hs=A*0;wavePERIOD=A*0;Uwave=A*0;waveANGLE=A*0;end
% 
% %SeaWaves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if computeSeaWave==1
% %hwave=ho;   %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% MASK=0*A+1;
% MASK(ho<=hwSea_lim | A==0 | VEG==1 | lev>=0)=0;
% [Uwave_sea,Tp_sea,Hsea,Fetch,kwave,PWsea]=SeaWaves(h,WINDdir,hwSea_lim,Trange,wind,MASK,64,dx);
% Uwave_sea=Uwave_sea.*(VEG==0); Hsea=Hsea.*(VEG==0); %vegetation effect
% [QsWslope_sea]=WaveSedimentTransport(Hsea,h,kwave,rhos,N,M,Tp_sea,dx,ss,ws1,hwSea_lim,fTide);
% QsWslope_sea(Hsea==0)=0;
% else;Uwave_sea=0*A;Tp_sea=0*A;Hsea=0*A;Fetch=0*A;QsWslope_sea=0*A;end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



