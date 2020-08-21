function [E,Ceq]=totalsedimenterosionSANDsine(h,hlimC,U,ss,d50,ws,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
g=9.81;
rho=1030;
rhos=2650;

manning=0.02;
Chezy=h.^(1/6)./manning;


% %Tides
% ncyc=10;
% Ceq_tide=0;
% for i=0:ncyc
% Ui=U*pi/2*sin(i/ncyc*pi/2);
% %Ui(Ui<0.2)=0;
% Ceq_tidei=0.05*Ui.^4./max(hlimC,h)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos; %kg/m/s%%
% %Ceq_tidei=0.05*Ui.^4./20./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos; %kg/m/s%
% Ceq_tide=Ceq_tide+1/(ncyc+1)*Ceq_tidei;
% end
% 

%Tides
ncyc=10;
Ceq_tide=0;
for i=0:ncyc
Ui=U*pi/2*sin(i/ncyc*pi/2);
%Ui(Ui<0.2)=0;
Ceq_tidei=Ui.^4; %kg/m/s
Ceq_tide=Ceq_tide+1/(ncyc+1)*Ceq_tidei;
end
Ceq_tide=0.05*Ceq_tide./max(hlimC,h)./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos;


%River
Ceq_river=0.05*UR.^4./h./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos*fMFriver.*fTide;


Ceq =Ceq_tide+Ceq_river;
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
% Chezy=h.^(1/8)./manning;
% Qs=0.05*(U*fUpeak).^5/fUpeak./(sqrt(g).*Chezy.^3*ss^2*d50)*rhos*Mfrequency; %kg/m/s

%meyepetermuller
% tau=1030*9.81.*h.^(-1/3).*manning.^2.*(U*fUpeak).^2;%tau=1030*0.0025*(U*fUpeak).^2;
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
% tau=1030*9.81.*h.^(-1/3).*manning.^2.*(Upeak).^2;%tau=1030*0.0025*(U*fUpeak).^2;
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
% tau=1030*9.81.*h.^(-1/3).*manning.^2.*(U*fUpeak).^2;
% Eu=me*max(0,tau-taucr)./taucr/fUpeak;
% qs=Eu;

%Di Silvio
% fUpeak=1;
% Ceq=0.00001*(fUpeak*U).^4/fUpeak./h*rhos;
% Eu=Ceq*(ws*3600*24); 
% qs=Eu;