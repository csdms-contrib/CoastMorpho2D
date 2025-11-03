function [h,ho,fTide,DHE,wl,wlo,hHR]=getwaterdepth(range,msl,z,kro,hpRIV,tempdeltaMSL);

% kro=0.01;
% hpRIV=[-3:0.01:5];
% z=hpRIV+0.0001;
% range=hpRIV*0+0.;
% tempdeltaMSL=0;
% msl=0;

%hpRIV is alwayd defined with respect to MSL and the tempMSL. It is becuase
%you awlays impose zero level at the open boundary and include msl and
%tempMLS in the water depth


ranged2=range/2;

%z=z-tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!
msl=msl+tempdeltaMSL;

%hR=max(0,msl+hpRIV-z);

maxHW=max(hpRIV,ranged2);
meanHW=hpRIV;
%minHW=hpRIV-range/2;% .*(hpRIV<range/2)
facTDRV=max(0,(ranged2-max(0,hpRIV-ranged2))./ranged2);


minHW=hpRIV-ranged2.*facTDRV;% .*(hpRIV<range/2)

rangeHW=min(range,maxHW-minHW);
%figure;imagesc(rangeHW);pause


%dHW=max(0,msl-z+hpRIV+range/2);%water depth at MHW
dHW=max(0,msl+maxHW-z);%water depth at MHW%April 24 2025

%fTide=min(1,max(10^-5,dHW./range));%hydroperiod
fTide=min(1,max(10^-5,dHW./rangeHW));%hydroperiod

%Option of using 2 or 3 water depth average
hw1=max(0,msl+maxHW-z);
hw2=max(0,msl+meanHW-z);
hw3=max(0,msl+minHW-z);
%h=hw2;

h=(hw1 +hw2 +hw3)/3; %for the flow
%h=(hw1*1 +hw3*9)/10; %for the flow
%h=(hw1 +hw2*3 +hw3)/5; %for the flow BESTA2025

%h=(hw1 +hw2*8 +hw3)/10; %for the flow
ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
h(h<kro)=kro;%reference water depth for flow calculation

%Efective water level disppalcement, creating the tidal prism
DHE=min(range,dHW);

%z=z+tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!
wl=z+h;
wlo=z+ho;



%hHR=(hw1+hw3)/2; %for the flow
%hHR=(hw1 +hw2 +hw3)/3; %for the flow %BESTA2025
%hHR=(hw1 +hw2*3 +hw3)/5; %for the flow BESTA2025


%hHR=(hw1 +hw2*2 +hw3)/4; %for the flow %besta sabaton 1 11 2025
%hHR=(hw1 +hw2*2 +hw3*2)/5; %for the flow
%hHR(hHR<kro)=kro;

hHR=NaN;



% figure;
% subplot(2,1,1);plot(hpRIV,minHW,hpRIV,maxHW,'.')
% subplot(2,1,2);plot(hpRIV,fTide)


%h3=(hw1 +hw2 +hw3)/3; %for the flow
%h3(h3<kro)=kro;%reference water depth for flow calculation
%hratio=max(1,h3./h);

% dBlo=-range/2+hpRIV;
% dBup=range/2+hpRIV;
% Zlev=z-msl;
% hratio=4*(Zlev-dBup).*(dBlo-Zlev)./(dBlo-dBup).^2;hratio(Zlev<dBlo)=0;hratio(Zlev>dBup)=0;%0.001; 
% hratio=hratio*(pi/2-1)+1;


% dBlo=-range/2+hpRIV;
% dBup=range/2+hpRIV;
% Zlev=z-msl;
% %hratio=z*0;
% hratio=(Zlev-dBlo)./range*2;
% hratio(Zlev>0)=1;
% hratio(Zlev<dBlo)=0;hratio(Zlev>dBup)=1;%0.001; 
% hratio=1+0.*hratio;
% %hratio=1+1*hratio;
% hratio(isnan(hratio))=1;
% 
% hratio=NaN;




% dBlo=-range/2+hpRIV;
% dBloi=dBlo-0.5;
% %dBup=range/2+hpRIV;
% Zlev=z-msl-dBlo;
% %hratio=z*0;
% hratio=(Zlev)/0.5+range+1;
% hratio(Zlev>dBlo)=1;
% hratio(Zlev<dBloi)=0;%hratio(Zlev>dBup)=1;%0.001; 
% hratio=1+1*hratio;

% figure
% imagesc(hratio)
% pause

% figure;
% plot(Zlev(:),hratio(:),'.')
% pause



%fTide=z*0+1;




% 
% figure
% subplot(3,1,1);imagesc(hw1');caxis([0 5.1])
% subplot(3,1,2);imagesc(hw2');caxis([0 5.1])
% subplot(3,1,3);imagesc(hw3');caxis([0 5.1])
% pause




%Original November 2019!!!!!!!!!!!!!
% function [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(range,Hsurge,msl,z,kro,hpRIV);
% DH=range+Hsurge; %total water displacement
% dHW=max(0,-z+hpRIV+msl+range/2+Hsurge);%water depth at MHW
% dtide=max(0,-z+hpRIV+msl+range/2);%water depth at MHW
% dsurge=max(0,-z+hpRIV+msl+Hsurge);%water depth at MHW
% fTide=min(1,max(10^-3,dHW/DH));%hydroperiod
% %h=0.5*(max(0,dHW)+max(0,dHW-DH)); %for the flow  THE ORGINIAL
% h=0.5*(max(0,dHW)+max(0,dHW-range)); %for the flow
% ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% h(h<kro)=kro;%reference water depth for flow calculation
% 
% % kl1=0.;
% % kl2=0.5;
% % a=find(h<kl2 & h>kl1);
% % h(a)=kl2/2+(h(a)-kl1)/(kl2-kl1)*kl2/2;
% % h(h<=kl1)=kl2/2;
% % 
% 
% %hi is the depth that gives the tidal prism %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% %ho=h;
% %ho(-z>(hpRIV+msl+range/2+Hsurge))=0; %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% 
% wl=z+h;
% wlo=z+ho;
% 
% % figure;
% % imagesc(ho)
% % pause