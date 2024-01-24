function [E,Ceq]=totalsedimenterosionSANDsine(h,hlimC,kro,MannS,rho,rhos,ss,d50,ws,fTide,computetidalcurrent,U,FcrUT,calculateshallowflow,US,hS,nSHALLOW,computeriver,UR,fMFriver);
g=9.81;
%nSHALLOW=1;
%hlimiterChezy=1;
%h=max(h,hlimiterChezy);%NEW February 2023
%Chezy=max(h,hlimC).^(1/6)./MannS;


Chezy=h.^(1/6)./MannS;
%Chezy=1./MannS;




%Tides%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetidalcurrent==1
UTmax=FcrUT*sqrt(9.8*h);
ncyc=10;
Ceq_tide=0;
for i=0:ncyc
Ui=U*pi/2*sin(i/ncyc*pi/2);
Ui=min(Ui,UTmax);
Ceq_tidei=0.05*Ui.^4./max(hlimC,h)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos; %kg/m/s
Ceq_tide=Ceq_tide+1/(ncyc+1)*Ceq_tidei;
end
else
Ceq_tide=0;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Shallow tidal flow%%%%%%%%%%%%%%%%
if calculateshallowflow==1
Chezy=hS.^(1/6)./MannS;
%Chezy=1./MannS;
Ceq_tideS=1/nSHALLOW *0.05*US.^4./max(hlimC,hS)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos;
%Ceq_tideS=1/nSHALLOW *0.05*US.^4./max(hlimC,h)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos;
else
Ceq_tideS=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%River%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeriver==1
FcrUR=0.5;
%hh=max(h,0.1);
UR=min(UR,FcrUR*sqrt(9.81*h));
UR=min(UR,3);
hlimC=0.01;
Chezy=max(1,h).^(1/6)./MannS;
%Chezy=5.^(1/6)./MannS;
%hlimitforU=1;
 % Ulimit=0.5;
% UR(h<hlimitforU & UR>Ulimit)=Ulimit+0.5*(UR(h<hlimitforU & UR>Ulimit)-Ulimit);   
%UR=min(UR,FcrUR*sqrt(9.81*h));
Ceq_river=0.05*UR.^4./max(hlimC,h)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos*fMFriver;%.*fTide;
else
Ceq_river=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ceq=Ceq_tide+Ceq_tideS+Ceq_river;
E=Ceq*(ws*3600*24);  































% Soulsby-van Rijn
% Ucr=0.19*d50^0.1*log10(4*max(h,0.1)/(2*d50));
% CD=(0.4./(log(h/0.006)-1)).^2;
% Ue=sqrt(Upeak.^2+0.018./CD.*Uwave_sea.^2);
% eps=max(0,(Ue-Ucr)/sqrt(ss*g*d50)).^2.4;
% Qsb=rhos*(eps*0.005.*h.*(d50./h).^1.2);
% Qss=rhos*(eps*0.012*d50*Ds.^-0.6);
% Qs_withoutU=fMFtide*(Qss+Qsb);

% Ue=Upeak;
% Qs_tide=QsSouslby(Ue,h,ss,g,d50,rhos,Ds);
% 
% CD=(0.4./(log(h/0.006)-1)).^2;
% Ue=0.018./CD.*Uwave_sea;
% Qs_wavesea=QsSouslby(Ue,h,ss,g,d50,rhos,Ds);
% 
% Qs_withoutU = fMFtide*Qs_tide + fMFsea*Qs_wavesea; 


















%figure;imagesc(max(0,(Upeak+0.4*Uwave_sea)));pause
%figure;imagesc(Ucr_sea);pause
%Eu(Utransport<=0)=0;

%Hengelung Hansen
% Chezy=h.^(1/8)./MannS;
% Qs=0.05*(U*fUpeak).^5/fUpeak./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos*Mfrequency; %kg/m/s

%meyepetermuller
% tau=1030*9.81.*h.^(-1/3).*MannS.^2.*(U*fUpeak).^2;%tau=1030*0.0025*(U*fUpeak).^2;
% theta=tau/((rhos-rho)*g*d50);
% Qs=sqrt(ss*g*d50^3)*8*max(0,(theta -0.06)).^1.5*rhos/fUpeak*Mfrequency; %0.047

%Soulsby van Rijn
%tcr=((rhos-rho)*g*d50)*0.06;
%Ucr=sqrt(tcr/(rho*0.0025));

% Ucr=0.19*d50^0.1*log10(4*h/(2*d50));Ucr(h<d50)=0;
% Ds=(g*ss/(10^-6)^2)^(1/3)*d50;
% Ass=(0.012*d50*Ds^-0.6)/(ss*g*d50)^1.2;
% Asb=(0.005*h.*(d50./h).^1.2)/(ss*g*d50)^1.2;
% As=Ass+Asb;
% Upeak=U*fUpeak;
% Urms=0;
% Ueff=max(0,sqrt((Upeak).^2 +Urms)-Ucr);
% Qs=As.*Upeak.*Ueff.^(2.4)*rhos/fUpeak*Mfrequency;

%Bijker
% Upeak=U*fUpeak;
% tau=1030*9.81.*h.^(-1/3).*MannS.^2.*(Upeak).^2;%tau=1030*0.0025*(U*fUpeak).^2;
% tauw=0;
% tautot=tau+tauw;
% Cb=2;
% muc=1;
% Qsb=Cb*d50*sqrt(tau*muc/rho).*exp(-0.27*(rhos-rho)*g.*d50./(muc.*tautot));
% Qss=1.83*Qsb*3;
% Qs=(Qsb+Qss)*rhos/fUpeak*Mfrequency;



% figure;imagesc(Qs');
% pause
















% taucr=0.3;
% me=50*10^-5*24*3600;
% tau=1030*9.81.*h.^(-1/3).*MannS.^2.*(U*fUpeak).^2;
% Eu=me*max(0,tau-taucr)./taucr/fUpeak;
% qs=Eu;

%Di Silvio
% fUpeak=1;
% Ceq=0.00001*(fUpeak*U).^4/fUpeak./h*rhos;
% Eu=Ceq*(ws*3600*24); 
% qs=Eu;