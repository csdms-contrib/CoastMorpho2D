function [h,ho,fTide,DHE,wl,wlo]=getwaterdepth(range,msl,z,kro,hpRIV,tempdeltaMSL);

z=z-tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!


dHW=max(0,msl-z+hpRIV+range/2);%water depth at MHW

fTide=min(1,max(10^-5,dHW./range));%hydroperiod

%Option of using 2 or 3 water depth average
hw1=max(0,dHW);
hw2=max(0,dHW-range/2);
hw3=max(0,dHW-range);
h=(hw1 +hw2 +hw3)/3; %for the flow
ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
h(h<kro)=kro;%reference water depth for flow calculation

%Efective water level disppalcement, creating the tidal prism
DHE=min(range,dHW);

wl=z+h;
wlo=z+ho;




% figure;
% imagesc(ho)
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