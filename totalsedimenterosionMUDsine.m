function  [E,Eflow,B]=totalsedimenterosionMUDsine(ho,h,taucr,taucrVEG,VEG,me,MannS,fTide,U,FcrUT,FcrUTveg,hlimC,calculateshallowflow,US,hS,nSHALLOW,computeSeaWaves,Uwave_seaI,Tp_seaI,computeSwellWave,Uwave,Tp_swell,computeRiver,fMFriver,UR,flowdestroyVEG,B);
          
%nSHALLOW=1;

taucro=U*0+taucr;
taucro(VEG==1)=taucrVEG;

%MannS=MannS*0+0.015;
Fcr=h*0;
Fcr(VEG==0)=FcrUT;
Fcr(VEG==1)=FcrUTveg;


Eflow=0;


UTmax=Fcr.*sqrt(9.8*ho);
%UTmax(h<=0.01 | VEG==1)=0;
%UTmax(h<=0.01 | VEG==1)=0;
ncyc=10;
Ei=0;
for i=0:ncyc
Utide=U*pi/2*sin(i/ncyc*pi/2);
Utide=min(Utide,UTmax);
%tauC=1030*9.81*MannS.^2.*max(h,hlimC).^(-1/3).*Utide.^2; 
tauC=1030*9.81*MannS.^2.*h.^(-1/3).*Utide.^2; 
%tauC=1030*9.81*MannS.^2.*Utide.^2; 
%tauC(h<0.5 & tauC>limitertauC)=limitertauC+(tauC(h<0.5 & tauC>limitertauC)-limitertauC)*0.1;%limiter for tay, since the smallest water depth kro = 5cm
Ei=Ei+1/(ncyc+1)*me.*(sqrt(1+(tauC./taucro).^2)-1);%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
%Ei=Ei+1/(ncyc+1)*me.*max(0,tauC-taucro)./taucro;%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
if flowdestroyVEG==1; B(tauC>taucro)=0; end
Eflow=Eflow+Ei;
end


%Shallow flow
if calculateshallowflow==1
tauC=1030*9.81*MannS.^2.*hS.^(-1/3).*US.^2; 
%tauC=1030*9.81*MannS.^2.*US.^2; 
%tauC(hS<0.5 & tauC>limitertauC)=limitertauC+(tauC(hS<0.5 & tauC>limitertauC)-limitertauC)*0.1;%limiter for tay, since the smallest water depth kro = 5cm
Ei=1/nSHALLOW*me.*(sqrt(1+(tauC./taucro).^2)-1);%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
%Ei=1/nSHALLOW*me.*max(0,tauC-taucro)./taucro;%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
Eflow=Eflow+Ei;
end


%Swell wave
if computeSwellWave==1
ko=1/1000;%
aw=Tp_swell.*Uwave/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;%fw=0.015;
tauWswell=0.5*1030*fw.*Uwave.^2;
%E=E+me*max(0,tauWswell-taucr)./taucro*fMFswell;
%E=E+me.*(sqrt(1+(tauWswell./taucro).^2)-1)l.*fTide; %OCIO AI TOLTO QUIESTI fTide IN APRIL 18th 2018
Eswellwave=me.*max(0,tauWswell-taucro)./taucro;
else
Eswellwave=0;
end


%Sea wave 
if computeSeaWaves==1
ko=1/1000;
Nh=size(Uwave_seaI,3);
Eseawave=0;
for i=1:Nh
Uwave_sea=max(0,Uwave_seaI(:,:,i));
Tp_sea=max(0,Tp_seaI(:,:,i));
aw=Tp_sea.*Uwave_sea/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;%fw=0.015;
tauWsea=0.5*1030*fw.*Uwave_sea.^2;
Ewtemp=me.*(sqrt(1+(tauWsea./taucro).^2)-1);
%Ewtemp=me.*max(0,tauWsea-taucro)./taucro;%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
Eseawave=Eseawave+Ewtemp;  %OCIO AI TOLTO QUIESTI IN APRIL 18th 2018
end
Eseawave=Eseawave/Nh;
else
Eseawave=0;
end


%River 
if computeRiver==1
%FcrUR=0.3;%UR=min(UR,FcrUR*sqrt(9.81*h));
%tauC=1030*9.81*MannS.^2.*max(h,hlimC).^(-1/3).*UR.^2;  
tauC=1030*9.81*MannS.^2.*h.^(-1/3).*UR.^2;  
%tauC=1030*9.81*MannS.^2.*UR.^2;  
%Ei=me.*max(0,tauC-taucro)./taucro*fMFriver.*fTide; 
Ei=me.*(sqrt(1+(tauC./taucro).^2)-1)*fMFriver;
Eflow=Eflow+Ei;
end




E=Eflow+Eseawave+Eswellwave;





















% tauCRiver=1030*0.04^2*9.81*h.^(-1/3).*UR.^2;
 

%E=E+1/(ncyc+1)*me./taucro.*(sqrt(taucr.^2+tauC.^2)-taucr);
%E=E+1/(ncyc+1)*me*max(0,tauC-taucr)./taucro;

%E=me./taucro.*(sqrt(taucr.^2+tauC.^2)-taucr)*fMFtide;


%NONCONTA!!!
%E=me*max(0,tauC-taucr)./taucro*fMFtide;% +me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%E=me*(sqrt(1+(tauC/taucr).^2)-1)*fMFtide;



%+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%%%
%E=h*NaN;
%Ec=me*max(0,tauC-taucr)./taucro*fMFtide;
%Ew=+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea;

%E=E*0;
%E(h<=0.1)=0;
% Ec(h<=0.1)=0;
% Ew(h<=0.1)=0;

%Ew=Ew*0;
%Ec=Ec*0;