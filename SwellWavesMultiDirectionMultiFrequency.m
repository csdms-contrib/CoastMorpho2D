function [H,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,d,ho,hwSwell_lim,Tperiodi,dx,periodic,angle,gridDIR,Ejonswap,nrefrac,wavediffraction);
%OCIO YOU IMPOSED PERIODIC==1 in chsoign the q+k later on (foe Delaware %trick)
nfrq=length(Tperiodi);

d(d<0 | A==0)=0;

%Energy balance of wind waves as a function of thebottom friction formulationR. Padilla-Hernandez) ´ , J. Monbaliu
%parameters
alpha=50; %adimensional
Cbr=0.73;
Cbed=0.038;
spread=15;
spread2=30;
spread3=45;
spread4=60;

Holateral=Ho;

    %DESCRIPTION
    %o is middle ray, 
    %i is left ray (angle turnerd clockwise, positive, increase the angle, whcih means will go toward lower left (more high angle).
    %j is right ray

%condtions on the boundaries (directinoal spding imposed at the b.c.)   
    facsprdo=1;
    facsprdi=0.;
    facsprdj=0.;
    facsprdi2=0.;
    facsprdj2=0.;
    facsprdi3=0;
    facsprdj3=0;
    facsprdi4=0;
    facsprdj4=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %   

anglei=NaN;anglej=NaN;anglei2=NaN;anglej2=NaN;anglei3=NaN;anglej3=NaN;anglei4=NaN;anglej4=NaN;
spdi=NaN;spdj=NaN;spdi2=NaN;spdj2=NaN;spdi3=NaN;spdj3=NaN;spdi4=NaN;spdj4=NaN;  

if nrefrac>=1
anglei=min(80,max(-80,angle+spread));spdi=anglei-angle;
anglej=min(80,max(-80,angle-spread));spdj=-(anglej-angle);
end
if nrefrac>=2
anglei2=min(80,max(-80,angle+spread2));spdi2=anglei2-angle-spdi;
anglej2=min(80,max(-80,angle-spread2));spdj2=-(anglej2-angle)-spdj;
end
if nrefrac>=3
anglei3=min(80,max(-80,angle+spread3));spdi3=anglei3-angle-spdi-spdi2;
anglej3=min(80,max(-80,angle-spread3));spdj3=-(anglej3-angle)-spdj-spdj2;
end
if nrefrac>=4
anglei4=min(80,max(-80,angle+spread4));spdi4=   anglei4-angle-spdi-spdi2-spdi3;
anglej4=min(80,max(-80,angle-spread4));spdj4=-(anglej4-angle)-spdj-spdj2-spdj3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Flip the domain upside down depending on the orientation of the grid   
if gridDIR==-1; 
    A=A(end:-1:1,:);
    AW=AW(end:-1:1,:);
    d=d(end:-1:1,:);
end



if nfrq>=1;kwave1=0*d+1;kwave1(d>0.2)=wavekgiulio(1/Tperiodi(1),d(d>0.2));end%kwave=wavek(1/Tp_swell,hwave);^Tperiodi(j)
if nfrq>=2;kwave2=0*d+1;kwave2(d>0.2)=wavekgiulio(1/Tperiodi(2),d(d>0.2));end%kwave=wavek(1/Tp_swell,hwave);^Tperiodi(j)
if nfrq>=3;kwave3=0*d+1;kwave3(d>0.2)=wavekgiulio(1/Tperiodi(3),d(d>0.2));end%kwave=wavek(1/Tp_swell,hwave);^Tperiodi(j)
if nfrq>=4;kwave4=0*d+1;kwave4(d>0.2)=wavekgiulio(1/Tperiodi(4),d(d>0.2));end%kwave=wavek(1/Tp_swell,hwave);^Tperiodi(j)
if nfrq>=5;kwave5=0*d+1;kwave5(d>0.2)=wavekgiulio(1/Tperiodi(5),d(d>0.2));end%kwave=wavek(1/Tp_swell,hwave);^Tperiodi(j)


Hcr=Cbr*d;Hcr(d<=0.1)=0;  
%IMPORTANT CHNAGE SEPT 2018!!!!!
Hcr(d<hwSwell_lim)=0;
    

if nfrq>=1;sigma1=2*pi/Tperiodi(1);end
if nfrq>=2;sigma2=2*pi/Tperiodi(2);end
if nfrq>=3;sigma3=2*pi/Tperiodi(3);end
if nfrq>=4;sigma4=2*pi/Tperiodi(4);end
if nfrq>=5;sigma5=2*pi/Tperiodi(5);end
if nfrq>=1;c1=(2*pi/Tperiodi(1))./kwave1;end
if nfrq>=2;c2=(2*pi/Tperiodi(2))./kwave2;end
if nfrq>=3;c3=(2*pi/Tperiodi(3))./kwave3;end
if nfrq>=4;c4=(2*pi/Tperiodi(4))./kwave4;end
if nfrq>=5;c5=(2*pi/Tperiodi(5))./kwave5;end
if nfrq>=1;cg1=0.5*c1.*(1+2*kwave1.*d./sinh(2*kwave1.*d));cg1(d<=0.1 | cg1<=0)=1; end
if nfrq>=2;cg2=0.5*c2.*(1+2*kwave2.*d./sinh(2*kwave2.*d));cg2(d<=0.1 | cg2<=0)=1; end 
if nfrq>=3;cg3=0.5*c3.*(1+2*kwave3.*d./sinh(2*kwave3.*d));cg3(d<=0.1 | cg3<=0)=1; end
if nfrq>=4;cg4=0.5*c4.*(1+2*kwave4.*d./sinh(2*kwave4.*d));cg4(d<=0.1 | cg4<=0)=1; end
if nfrq>=5;cg5=0.5*c5.*(1+2*kwave5.*d./sinh(2*kwave5.*d));cg5(d<=0.1 | cg5<=0)=1; end
if nfrq>=1;bedfr1=1*(sigma1./(9.81*sinh(kwave1.*d))).^2;bedfr1(d<=0.1)=1;end;% there use to be a Cbed coefficinet (now 1). substided in waveEVOLUTIONstep.m
if nfrq>=2;bedfr2=1*(sigma2./(9.81*sinh(kwave2.*d))).^2;bedfr2(d<=0.1)=1;end
if nfrq>=3;bedfr3=1*(sigma3./(9.81*sinh(kwave3.*d))).^2;bedfr3(d<=0.1)=1;end
if nfrq>=4;bedfr4=1*(sigma4./(9.81*sinh(kwave4.*d))).^2;bedfr4(d<=0.1)=1;end
if nfrq>=5;bedfr5=1*(sigma5./(9.81*sinh(kwave5.*d))).^2;bedfr5(d<=0.1)=1;end


%Refraction coefficients. There is the extra 2 for calibration
if nrefrac>=1
    if nfrq>=1;dCy1=([c1(:,2:end) c1(:,end)]-[c1(:,1) c1(:,1:end-1)])/2./c1.*(d>1)*1;dCx1=([c1(2:end,:); c1(end,:)]-[c1(1,:); c1(1:end-1,:)])/2./c1.*(d>1)*1;end
    if nfrq>=2;dCy2=([c2(:,2:end) c2(:,end)]-[c2(:,1) c2(:,1:end-1)])/2./c2.*(d>1)*1;dCx2=([c2(2:end,:); c2(end,:)]-[c2(1,:); c2(1:end-1,:)])/2./c2.*(d>1)*1;end
    if nfrq>=3;dCy3=([c3(:,2:end) c3(:,end)]-[c3(:,1) c3(:,1:end-1)])/2./c3.*(d>1)*1;dCx3=([c3(2:end,:); c3(end,:)]-[c3(1,:); c3(1:end-1,:)])/2./c3.*(d>1)*1;end
    if nfrq>=4;dCy4=([c4(:,2:end) c4(:,end)]-[c4(:,1) c4(:,1:end-1)])/2./c4.*(d>1)*1;dCx4=([c4(2:end,:); c4(end,:)]-[c4(1,:); c4(1:end-1,:)])/2./c4.*(d>1)*1;end
    if nfrq>=5;dCy5=([c5(:,2:end) c5(:,end)]-[c5(:,1) c5(:,1:end-1)])/2./c5.*(d>1)*1;dCx5=([c5(2:end,:); c5(end,:)]-[c5(1,:); c5(1:end-1,:)])/2./c5.*(d>1)*1;end
end


H=zeros(N,M);H(end,:)=Ho;
waveANGLE=d*NaN;%just to initialize
wavePERIOD=d*0;%just to initialize
PWswell=d*0;   


p = [1:M]';%exclude the NOLAND CELLS (A==0)
%%%%%%%%%%%%%%%%%%INITIALIZE THE INDEXES AND VECTOES
ii=[];jj=[];
for k=[-1 1];   
%if periodic==0;     q=p+k;if k==1;a=[1:M-1];elseif k==-1;a=[2:M];end 
%elseif periodic==1; q=1+mod(p+k-1,M);a=[1:M];
%end
q=1+mod(p+k-1,M);a=[1:M];%%KINDA force the periodic option in chosiing the cell (Delaware wave diffusion)

ii=[ii;p(a)]; jj=[jj;q(a)];%gain from the neigborh cell
end
ii=[ii;p];jj=[jj;p];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ALLfacsprd=[facsprdo,facsprdi,facsprdj,facsprdi2,facsprdj2,facsprdi3,facsprdj3,facsprdi4,facsprdj4];
    ALLangle=[angle,anglei,anglej,anglei2,anglej2,anglei3,anglej3,anglei4,anglej4];
    ALLspd=[spdi,spdj,spdi2,spdj2,spdi3,spdj3,spdi4,spdj4];  
    

%The initial E value at the sea boundary
if nfrq>=1;E1=H(end,:).^2/16*Ejonswap(1);end
if nfrq>=2;E2=H(end,:).^2/16*Ejonswap(2);end
if nfrq>=3;E3=H(end,:).^2/16*Ejonswap(3);end
if nfrq>=4;E4=H(end,:).^2/16*Ejonswap(4);end
if nfrq>=5;E5=H(end,:).^2/16*Ejonswap(5);end

if nfrq==1;E=E1;end
if nfrq==2;E=E1+E2;end
if nfrq==3;E=E1+E2+E3;end
if nfrq==4;E=E1+E2+E3+E4;end
if nfrq==5;E=E1+E2+E3+E4+E5;end

    
%Partition the wave heigth in the boundary condtions    

if nfrq>=1;
    Eo=E1*facsprdo;
    if nrefrac>=1;    Ei=E1*facsprdi;      Ej=E1*facsprdj;    end
    if nrefrac>=2;    Ei2=E1*facsprdi2;    Ej2=E1*facsprdj2;    end
    if nrefrac>=3;    Ei3=E1*facsprdi3;    Ej3=E1*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E1*facsprdi4;    Ej4=E1*facsprdj4;    end
if nrefrac==0;ALLE1=Eo;elseif nrefrac==1;ALLE1=[Eo;Ei;Ej];elseif nrefrac==2;ALLE1=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE1=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE1=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end
end

if nfrq>=2;
    Eo=E2*facsprdo;
    if nrefrac>=1;    Ei=E2*facsprdi;      Ej=E2*facsprdj;    end
    if nrefrac>=2;    Ei2=E2*facsprdi2;    Ej2=E2*facsprdj2;    end
    if nrefrac>=3;    Ei3=E2*facsprdi3;    Ej3=E2*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E2*facsprdi4;    Ej4=E2*facsprdj4;    end
if nrefrac==0;ALLE2=Eo;elseif nrefrac==1;ALLE2=[Eo;Ei;Ej];elseif nrefrac==2;ALLE2=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE2=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE2=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end
end

if nfrq>=3;
    Eo=E3*facsprdo;
    if nrefrac>=1;    Ei=E3*facsprdi;      Ej=E3*facsprdj;    end
    if nrefrac>=2;    Ei2=E3*facsprdi2;    Ej2=E3*facsprdj2;    end
    if nrefrac>=3;    Ei3=E3*facsprdi3;    Ej3=E3*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E3*facsprdi4;    Ej4=E3*facsprdj4;    end
if nrefrac==0;ALLE3=Eo;elseif nrefrac==1;ALLE3=[Eo;Ei;Ej];elseif nrefrac==2;ALLE3=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE3=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE3=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end
end

if nfrq>=4;
    Eo=E4*facsprdo;
    if nrefrac>=1;    Ei=E4*facsprdi;      Ej=E4*facsprdj;    end
    if nrefrac>=2;    Ei2=E4*facsprdi2;    Ej2=E4*facsprdj2;    end
    if nrefrac>=3;    Ei3=E4*facsprdi3;    Ej3=E4*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E4*facsprdi4;    Ej4=E4*facsprdj4;    end
if nrefrac==0;ALLE4=Eo;elseif nrefrac==1;ALLE4=[Eo;Ei;Ej];elseif nrefrac==2;ALLE4=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE4=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE4=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end
end

if nfrq>=5;
    Eo=E5*facsprdo;
    if nrefrac>=1;    Ei=E5*facsprdi;      Ej=E5*facsprdj;    end
    if nrefrac>=2;    Ei2=E5*facsprdi2;    Ej2=E5*facsprdj2;    end
    if nrefrac>=3;    Ei3=E5*facsprdi3;    Ej3=E5*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E5*facsprdi4;    Ej4=E5*facsprdj4;    end
if nrefrac==0;ALLE5=Eo;elseif nrefrac==1;ALLE5=[Eo;Ei;Ej];elseif nrefrac==2;ALLE5=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE5=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE5=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACTUAL LOOP from the wave boundary and up following the waves  
for i=N-1:-1:1  
    
if nfrq==1;sigmatot=2*pi*(max(0,E1./E)/Tperiodi(1));sigmatot(E==0)=0;end%   *(E1/Tperiodi(1)+E2/Tperiodi(2))./max(0.00001,E);
if nfrq==2;sigmatot=2*pi*(max(0,E1./E)/Tperiodi(1)+max(0,E2./E)/Tperiodi(2));sigmatot(E==0)=0;end%   *(E1/Tperiodi(1)+E2/Tperiodi(2))./max(0.00001,E);
if nfrq==3;sigmatot=2*pi*(max(0,E1./E)/Tperiodi(1)+max(0,E2./E)/Tperiodi(2)+max(0,E3./E)/Tperiodi(3));sigmatot(E==0)=0;end%   *(E1/Tperiodi(1)+E2/Tperiodi(2))./max(0.00001,E);
if nfrq==4;sigmatot=2*pi*(max(0,E1./E)/Tperiodi(1)+max(0,E2./E)/Tperiodi(2)+max(0,E3./E)/Tperiodi(3)+max(0,E4./E)/Tperiodi(4));sigmatot(E==0)=0;end%   *(E1/Tperiodi(1)+E2/Tperiodi(2))./max(0.00001,E);
if nfrq==5;sigmatot=2*pi*(max(0,E1./E)/Tperiodi(1)+max(0,E2./E)/Tperiodi(2)+max(0,E3./E)/Tperiodi(3)+max(0,E4./E)/Tperiodi(4)+max(0,E5./E)/Tperiodi(5));sigmatot(E==0)=0;end%   *(E1/Tperiodi(1)+E2/Tperiodi(2))./max(0.00001,E);

wavePERIOD(i,sigmatot>0)=2*pi./sigmatot(sigmatot>0);


if nfrq==1;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1;end
if nfrq==2;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2;end
if nfrq==3;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2+cg3(i,:)*1030*9.8.*E3;end
if nfrq==4;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2+cg3(i,:)*1030*9.8.*E3+cg4(i,:)*1030*9.8.*E4;end
if nfrq==5;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2+cg3(i,:)*1030*9.8.*E3+cg4(i,:)*1030*9.8.*E4+cg5(i,:)*1030*9.8.*E5;end

% if nfrq==1;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1;end
% if nfrq==2;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2;end
% if nfrq==3;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2+cg3(i,:)*1030*9.8.*E3;end
% if nfrq==4;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1/Tperiodi(1)+cg2(i,:)*1030*9.8.*E2/Tperiodi(2)+cg3(i,:)*1030*9.8.*E3/Tperiodi(3)+cg4(i,:)*1030*9.8.*E4/Tperiodi(4);end
% if nfrq==5;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1/Tperiodi(1)+cg2(i,:)*1030*9.8.*E2/Tperiodi(2)+cg3(i,:)*1030*9.8.*E3/Tperiodi(3)+cg4(i,:)*1030*9.8.*E4/Tperiodi(4)+cg5(i,:)*1030*9.8.*E5/Tperiodi(5);end

% Tcutoff=4;
% %if nfrq==1;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1;end
% %if nfrq==2;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2;end
% %if nfrq==3;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1+cg2(i,:)*1030*9.8.*E2+cg3(i,:)*1030*9.8.*E3;end
% if nfrq==4;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1*(Tperiodi(1)>Tcutoff)+cg2(i,:)*1030*9.8.*E2*(Tperiodi(2)>Tcutoff)+cg3(i,:)*1030*9.8.*E3*(Tperiodi(3)>Tcutoff)+cg4(i,:)*1030*9.8.*E4*(Tperiodi(4)>Tcutoff);end
% if nfrq==5;PWswell(i,:)=cg1(i,:)*1030*9.8.*E1*(Tperiodi(1)>Tcutoff)+cg2(i,:)*1030*9.8.*E2*(Tperiodi(2)>Tcutoff)+cg3(i,:)*1030*9.8.*E3*(Tperiodi(3)>Tcutoff)+cg4(i,:)*1030*9.8.*E4*(Tperiodi(4)>Tcutoff)+cg5(i,:)*1030*9.8.*E5*(Tperiodi(5)>Tcutoff);end


% if nfrq>=1;[E1,ALLE1]=waveEVOLUTIONstep(E,E1,ALLE1,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i+1,:),kwave1(i+1,:),c1,cg1(i,:),cg1(i+1,:),bedfr1(i+1,:),sigma1,Tperiodi(1),dx,periodic,Holateral*sqrt(Ejonswap(1)),sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx1(i+1,:),dCy1(i+1,:));end
% 
% if nfrq>=2;[E2,ALLE2]=waveEVOLUTIONstep(E,E2,ALLE2,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i+1,:),kwave2(i+1,:),c2,cg2(i,:),cg2(i+1,:),bedfr2(i+1,:),sigma2,Tperiodi(2),dx,periodic,Holateral*sqrt(Ejonswap(2)),sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx2(i+1,:),dCy2(i+1,:)); end
% 
% if nfrq>=3;   [E3,ALLE3]=waveEVOLUTIONstep(E,E3,ALLE3,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i+1,:),kwave3(i+1,:),c3,cg3(i,:),cg3(i+1,:),bedfr3(i+1,:),sigma3,Tperiodi(3),dx,periodic,Holateral*sqrt(Ejonswap(3)),sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx3(i+1,:),dCy3(i+1,:));  end
%  
% if nfrq>=4;  [E4,ALLE4]=waveEVOLUTIONstep(E,E4,ALLE4,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i+1,:),kwave4(i+1,:),c4,cg4(i,:),cg4(i+1,:),bedfr4(i+1,:),sigma4,Tperiodi(4),dx,periodic,Holateral*sqrt(Ejonswap(4)),sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx4(i+1,:),dCy4(i+1,:));  end
% 
% if nfrq>=5;  [E5,ALLE5]=waveEVOLUTIONstep(E,E5,ALLE5,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i+1,:),kwave5(i+1,:),c5,cg5(i,:),cg5(i+1,:),bedfr5(i+1,:),sigma5,Tperiodi(5),dx,periodic,Holateral*sqrt(Ejonswap(5)),sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx4(i+1,:),dCy4(i+1,:));  end


if nfrq>=1;[E1,ALLE1,dPW]=waveEVOLUTIONstep(E,E1,ALLE1,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave1(i,:),c1,cg1(i,:),cg1(i+1,:),bedfr1(i,:),sigma1,Tperiodi(1),dx,periodic,Holateral*sqrt(Ejonswap(1)),sigmatot,...
       ALLangle,ALLfacsprd,ALLspd,dCx1(i+1,:),dCy1(i+1,:));end

if nfrq>=2;[E2,ALLE2,dPW]=waveEVOLUTIONstep(E,E2,ALLE2,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave2(i,:),c2,cg2(i,:),cg2(i+1,:),bedfr2(i,:),sigma2,Tperiodi(2),dx,periodic,Holateral*sqrt(Ejonswap(2)),sigmatot,...
       ALLangle,ALLfacsprd,ALLspd,dCx2(i+1,:),dCy2(i+1,:)); end

if nfrq>=3;[E3,ALLE3,dPW]=waveEVOLUTIONstep(E,E3,ALLE3,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave3(i,:),c3,cg3(i,:),cg3(i+1,:),bedfr3(i,:),sigma3,Tperiodi(3),dx,periodic,Holateral*sqrt(Ejonswap(3)),sigmatot,...
       ALLangle,ALLfacsprd,ALLspd,dCx3(i+1,:),dCy3(i+1,:));  end
 
if nfrq>=4;[E4,ALLE4,dPW]=waveEVOLUTIONstep(E,E4,ALLE4,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave4(i,:),c4,cg4(i,:),cg4(i+1,:),bedfr4(i,:),sigma4,Tperiodi(4),dx,periodic,Holateral*sqrt(Ejonswap(4)),sigmatot,...
       ALLangle,ALLfacsprd,ALLspd,dCx4(i+1,:),dCy4(i+1,:));  end

if nfrq>=5;[E5,ALLE5,dPW]=waveEVOLUTIONstep(E,E5,ALLE5,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave5(i,:),c5,cg5(i,:),cg5(i+1,:),bedfr5(i,:),sigma5,Tperiodi(5),dx,periodic,Holateral*sqrt(Ejonswap(5)),sigmatot,...
       ALLangle,ALLfacsprd,ALLspd,dCx4(i+1,:),dCy4(i+1,:));  end

if nfrq==1;E=E1;end
if nfrq==2;E=E1+E2;end
if nfrq==3;E=E1+E2+E3;end
if nfrq==4;E=E1+E2+E3+E4;end
if nfrq==5;E=E1+E2+E3+E4+E5;end


    if nfrq>=1; ffE1=max(0,E1./E);end
    if nfrq>=2; ffE2=max(0,E2./E);end
    if nfrq>=3; ffE3=max(0,E3./E);end
    if nfrq>=4; ffE4=max(0,E4./E);end
    if nfrq>=5; ffE5=max(0,E5./E);end
    
    

    if nfrq>=1;
     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f1Eo=max(0,ALLE1(1,:)./E1);     f1Ei=max(0,ALLE1(2,:)./E1);          f1Ej=max(0,ALLE1(3,:)./E1); 
        if nrefrac>=2;        f1Ei2=max(0,ALLE1(4,:)./E1);        f1Ej2=max(0,ALLE1(5,:)./E1);         end
        if nrefrac>=3;        f1Ei3=max(0,ALLE1(6,:)./E1);        f1Ej3=max(0,ALLE1(7,:)./E1);         end
        if nrefrac>=4;        f1Ei4=max(0,ALLE1(8,:)./E1);        f1Ej4=max(0,ALLE1(9,:)./E1);         end
     end
    end  
    if nfrq>=2;
     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f2Eo=max(0,ALLE2(1,:)./E2);     f2Ei=max(0,ALLE2(2,:)./E2);          f2Ej=max(0,ALLE2(3,:)./E2); 
        if nrefrac>=2;        f2Ei2=max(0,ALLE2(4,:)./E2);        f2Ej2=max(0,ALLE2(5,:)./E2);         end
        if nrefrac>=3;        f2Ei3=max(0,ALLE2(6,:)./E2);        f2Ej3=max(0,ALLE2(7,:)./E2);         end
        if nrefrac>=4;        f2Ei4=max(0,ALLE2(8,:)./E2);        f2Ej4=max(0,ALLE2(9,:)./E2);         end
     end
    end   
    if nfrq>=3;
     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f3Eo=max(0,ALLE3(1,:)./E3);     f3Ei=max(0,ALLE3(2,:)./E3);          f3Ej=max(0,ALLE3(3,:)./E3); 
        if nrefrac>=2;        f3Ei2=max(0,ALLE3(4,:)./E3);        f3Ej2=max(0,ALLE3(5,:)./E3);         end
        if nrefrac>=3;        f3Ei3=max(0,ALLE3(6,:)./E3);        f3Ej3=max(0,ALLE3(7,:)./E3);         end
        if nrefrac>=4;        f3Ei4=max(0,ALLE3(8,:)./E3);        f3Ej4=max(0,ALLE3(9,:)./E3);         end
     end
    end  
    if nfrq>=4;
     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f4Eo=max(0,ALLE4(1,:)./E4);     f4Ei=max(0,ALLE4(2,:)./E4);          f4Ej=max(0,ALLE4(3,:)./E4); 
        if nrefrac>=2;        f4Ei2=max(0,ALLE4(4,:)./E4);        f4Ej2=max(0,ALLE4(5,:)./E4);         end
        if nrefrac>=3;        f4Ei3=max(0,ALLE4(6,:)./E4);        f4Ej3=max(0,ALLE4(7,:)./E4);         end
        if nrefrac>=4;        f4Ei4=max(0,ALLE4(8,:)./E4);        f4Ej4=max(0,ALLE4(9,:)./E4);         end
     end
    end  
    if nfrq>=5;
     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f5Eo=max(0,ALLE5(1,:)./E5);     f5Ei=max(0,ALLE5(2,:)./E5);          f5Ej=max(0,ALLE5(3,:)./E5); 
        if nrefrac>=2;        f5Ei2=max(0,ALLE5(4,:)./E5);        f5Ej2=max(0,ALLE5(5,:)./E5);         end
        if nrefrac>=3;        f5Ei3=max(0,ALLE5(6,:)./E5);        f5Ej3=max(0,ALLE5(7,:)./E5);         end
        if nrefrac>=4;        f5Ei4=max(0,ALLE5(8,:)./E5);        f5Ej4=max(0,ALLE5(9,:)./E5);         end
     end
    end
    
   
    
   wANGtot=0;        
    if nfrq>=1;
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f1Eo+anglei*f1Ei+anglej*f1Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f1Ei2+anglej2*f1Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f1Ei3+anglej3*f1Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f1Ei4+anglej4*f1Ej4;        end
        wANGtot=wANGtot+wANG.*ffE1;
    end 
    if nfrq>=2;    
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f2Eo+anglei*f2Ei+anglej*f2Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f2Ei2+anglej2*f2Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f2Ei3+anglej3*f2Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f2Ei4+anglej4*f2Ej4;        end
        wANGtot=wANGtot+wANG.*ffE2;
    end    
    if nfrq>=3;    
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f3Eo+anglei*f3Ei+anglej*f3Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f3Ei2+anglej2*f3Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f3Ei3+anglej3*f3Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f3Ei4+anglej4*f3Ej4;        end
        wANGtot=wANGtot+wANG.*ffE3;
    end      
    if nfrq>=4; 
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f4Eo+anglei*f4Ei+anglej*f4Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f4Ei2+anglej2*f4Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f4Ei3+anglej3*f4Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f4Ei4+anglej4*f4Ej4;        end
        wANGtot=wANGtot+wANG.*ffE4;
    end      
    if nfrq>=5;
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f5Eo+anglei*f5Ei+anglej*f5Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f5Ei2+anglej2*f5Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f5Ei3+anglej3*f5Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f5Ei4+anglej4*f5Ej4;        end
        wANGtot=wANGtot+wANG.*ffE5;
    end
        
        
        waveANGLE(i,:)=wANGtot;  
        
    %Diffract/Diffuse
    if wavediffraction==1; 
        %     if nfrq==1;Mc=c1(i+1,:).*ffE1;Mcg=cg1(i+1,:).*ffE1;end
        %     if nfrq==2;Mc=c1(i+1,:).*ffE1+c2(i+1,:).*ffE2;Mcg=cg1(i+1,:).*ffE1+cg2(i+1,:).*ffE2;end
        %     if nfrq==3;Mc=c1(i+1,:).*ffE1+c2(i+1,:).*ffE2+c3(i+1,:).*ffE3;Mcg=cg1(i+1,:).*ffE1+cg2(i+1,:).*ffE2+cg3(i+1,:).*ffE3;end
        %     if nfrq==4;Mc=c1(i+1,:).*ffE1+c2(i+1,:).*ffE2+c3(i+1,:).*ffE3+c4(i+1,:).*ffE4;Mcg=cg1(i+1,:).*ffE1+cg2(i+1,:).*ffE2+cg3(i+1,:).*ffE3+cg4(i+1,:).*ffE4;end
        %     if nfrq==5;Mc=c1(i+1,:).*ffE1+c2(i+1,:).*ffE2+c3(i+1,:).*ffE3+c4(i+1,:).*ffE4+c5(i+1,:).*ffE5;Mcg=cg1(i+1,:).*ffE1+cg2(i+1,:).*ffE2+cg3(i+1,:).*ffE3+cg4(i+1,:).*ffE4+cg5(i+1,:).*ffE5;end     
    Mc=c1(i,:);Mcg=cg1(i,:);
    [E]=wavediffusionstep(E,i,Mc,Mcg,d(i,:),ho(i,:),sigmatot,dx,angle,Cbr,Cbed,alpha,angle,N,M,A,AW(i,:),AW(i+1,:),p,ii,jj,periodic,Holateral); 
    end
    
    %Recalculate the wave height and set breaking criterion check
    H(i,:)=4*sqrt(max(0,E));
    a=find(H(i,:)>Hcr(i,:));
    H(i,a)=Hcr(i,H(i,:)>Hcr(i,:));
    E(a)=Hcr(i,a).^2/16; 
    
    %breaking for longshore transport
    a=find(H(i+1,:)>Hcr(i,:));%waves will break!
    deltaPW(i+1,a)=max(0,H(i+1,a).^2.*cg1(i+1,a)-Hcr(i,a).^2.*cg1(i,a))/dx*1030*9.8/16;
    
    
  if nfrq>=1;
    E1=ffE1.*E;
    ALLE1(1,:)=E1.*f1Eo;
    if nrefrac>=1;        ALLE1(2,:)=E1.*f1Ei;         ALLE1(3,:)=E1.*f1Ej;        end
    if nrefrac>=2;        ALLE1(4,:)=E1.*f1Ei2;        ALLE1(5,:)=E1.*f1Ej2;       end
    if nrefrac>=3;        ALLE1(6,:)=E1.*f1Ei3;        ALLE1(7,:)=E1.*f1Ej3;       end
    if nrefrac>=4;        ALLE1(8,:)=E1.*f1Ei4;        ALLE1(9,:)=E1.*f1Ej4;       end
  end
  if nfrq>=2;  
    E2=ffE2.*E;
    ALLE2(1,:)=E2.*f2Eo;
    if nrefrac>=1;        ALLE2(2,:)=E2.*f2Ei;         ALLE2(3,:)=E2.*f2Ej;        end
    if nrefrac>=2;        ALLE2(4,:)=E2.*f2Ei2;        ALLE2(5,:)=E2.*f2Ej2;       end
    if nrefrac>=3;        ALLE2(6,:)=E2.*f2Ei3;        ALLE2(7,:)=E2.*f2Ej3;       end
    if nrefrac>=4;        ALLE2(8,:)=E2.*f2Ei4;        ALLE2(9,:)=E2.*f2Ej4;       end
  end
  if nfrq>=3;   
    E3=ffE3.*E;
    ALLE3(1,:)=E3.*f3Eo;
    if nrefrac>=1;        ALLE3(2,:)=E3.*f3Ei;         ALLE3(3,:)=E3.*f3Ej;        end
    if nrefrac>=2;        ALLE3(4,:)=E3.*f3Ei2;        ALLE3(5,:)=E3.*f3Ej2;       end
    if nrefrac>=3;        ALLE3(6,:)=E3.*f3Ei3;        ALLE3(7,:)=E3.*f3Ej3;       end
    if nrefrac>=4;        ALLE3(8,:)=E3.*f3Ei4;        ALLE3(9,:)=E3.*f3Ej4;       end
  end
  if nfrq>=4;
    E4=ffE4.*E; 
    ALLE4(1,:)=E4.*f4Eo;
    if nrefrac>=1;        ALLE4(2,:)=E4.*f4Ei;         ALLE4(3,:)=E4.*f4Ej;        end
    if nrefrac>=2;        ALLE4(4,:)=E4.*f4Ei2;        ALLE4(5,:)=E4.*f4Ej2;       end
    if nrefrac>=3;        ALLE4(6,:)=E4.*f4Ei3;        ALLE4(7,:)=E4.*f4Ej3;       end
    if nrefrac>=4;        ALLE4(8,:)=E4.*f4Ei4;        ALLE4(9,:)=E4.*f4Ej4;       end
  end
  if nfrq>=5;
    E5=ffE5.*E;
    ALLE5(1,:)=E5.*f5Eo;
    if nrefrac>=1;        ALLE5(2,:)=E5.*f5Ei;         ALLE5(3,:)=E5.*f5Ej;        end
    if nrefrac>=2;        ALLE5(4,:)=E5.*f5Ei2;        ALLE5(5,:)=E5.*f5Ej2;       end
    if nrefrac>=3;        ALLE5(6,:)=E5.*f5Ei3;        ALLE5(7,:)=E5.*f5Ej3;       end
    if nrefrac>=4;        ALLE5(8,:)=E5.*f5Ei4;        ALLE5(9,:)=E5.*f5Ej4;       end
  end
       % waveANGLE(N-i,:)=wANG;    
end;%end of wave evolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uwave=(pi*H./wavePERIOD./sinh(kwave1.*d));Uwave(d<=0)=0;
%H=sqrt(H1.^2+H2.^2);
%wavePERIOD=1./((H1.^2/Tperiodi(1)+H2.^2/Tperiodi(2))./H.^2);

kwave=kwave1;

%flip the domain back
if gridDIR==-1; 
    H=H(end:-1:1,:);
    wavePERIOD=wavePERIOD(end:-1:1,:);
    waveANGLE=waveANGLE(end:-1:1,:);
    PWswell=PWswell(end:-1:1,:);
    kwave=kwave(end:-1:1,:);
    Uwave=Uwave(end:-1:1,:);
    deltaPW=deltaPW(end:-1:1,:);
end







    %repartition the energy
%     Eo=E;
%     if nrefrac>=1;
%         Eo=E.*fEo;        Ei=E.*fEi;        Ej=E.*fEj;
%         if nrefrac>=2;        Ei2=E.*fEi2;        Ej2=E.*fEj2;          end
%         if nrefrac>=3;        Ei3=E.*fEi3;        Ej3=E.*fEj3;            end
%         if nrefrac>=4;        Ei4=E.*fEi4;        Ej4=E.*fEj4;            end
%     end
%     
% if nrefrac==0;ALLE=Eo;elseif nrefrac==1;ALLE=[Eo;Ei;Ej];elseif nrefrac==2;ALLE=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
% elseif nrefrac==4;ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end 






% if nrefrac>=1
% %dCy=([c(:,2:end) c(:,end)]-[c(:,1) c(:,1:end-1)])/2./c.*(d>0.1);
% %dCx=([c(2:end,:); c(end,:)]-[c(1,:); c(1:end-1,:)])/2./c.*(d>0.1);
% dCi=sin((angle0+anglei0)/2/180*pi)*dCx-cos((angle0+anglei0)/2/180*pi)*dCy;
% dCj=sin((angle0+anglej0)/2/180*pi)*dCx-cos((angle0+anglej0)/2/180*pi)*dCy;
% end
% if nrefrac>=2
% dCi2=sin((anglei0+anglei20)/2/180*pi)*dCx-cos((anglei0+anglei20)/2/180*pi)*dCy;
% dCj2=sin((anglej0+anglej20)/2/180*pi)*dCx-cos((anglej0+anglej20)/2/180*pi)*dCy;
% end
% if nrefrac>=3
% dCi3=sin((anglei20+anglei30)/2/180*pi)*dCx-cos((anglei20+anglei30)/2/180*pi)*dCy;
% dCj3=sin((anglej20+anglej30)/2/180*pi)*dCx-cos((anglej20+anglej30)/2/180*pi)*dCy;
% end
% if nrefrac>=4
% dCi4=sin((anglei30+anglei40)/2/180*pi)*dCx-cos((anglei30+anglei40)/2/180*pi)*dCy;
% dCj4=sin((anglej30+anglej40)/2/180*pi)*dCx-cos((anglej30+anglej40)/2/180*pi)*dCy;
% end




% %ACTUAL LOOP
% for i=N-1:-1:1
%     
%     
%     Eo=E*(1-1*facspread);
%     %Ei=E*facspread;
%     Ej=E*facspread;
%        
%     s=floor(tgangle*(i+1));Eo=circshift(Eo,[0 s*angledir]);
%     [Eo]=wavepropagationstep(Eo,i,c,cg,d,kwave,sigma,T,dx,angle,Cbr,Cbed,alpha,angle0,N,M,A,AW,p,ii,jj,periodic,Holateral*sqrt(1-1*facspread)); 
%     s=floor(tgangle*i); Eo=circshift(Eo,[0 -s*angledir]);
%     H(i,:)=4*sqrt(max(0,Eo));a=find(H(i,:)>Hcr(i,:));H(i,a)=Hcr(i,a);Eo(a)=Hcr(i,a).^2/16;
% %     
% %     s=floor(tganglei*(i+1));Ei=circshift(Ei,[0 s*anglediri]);
% %     [Ei]=wavepropagationstep(Ei,i,ci,cgi,di,kwavei,sigma,T,dx,anglei,Cbr,Cbed,alpha,anglei0,N,M,A,AWi,p,ii,jj,periodic,Holateral*sqrt(facspread)); 
% %     s=floor(tganglei*i); Ei=circshift(Ei,[0 -s*anglediri]);
%     
%     s=floor(tganglej*(i+1));Ej=circshift(Ej,[0 s*angledirj]);
%     [Ej]=wavepropagationstep(Ej,i,cj,cgj,dj,kwavej,sigma,T,dx,anglej,Cbr,Cbed,alpha,anglej0,N,M,A,AWj,p,ii,jj,periodic,Holateral*sqrt(facspread)); 
%     s=floor(tganglej*i); Ej=circshift(Ej,[0 -s*angledirj]);
%     H(i,:)=4*sqrt(max(0,Ej));a=find(H(i,:)>Hcr(i,:));H(i,a)=Hcr(i,a);Ej(a)=Hcr(i,a).^2/16;
% %     
%     %E=Ei+Ej;
%     E=Eo+Ej;   
%     
%     
%     %recaluclate the wave height
%     H(i,:)=4*sqrt(max(0,E));
% % 
% %         %breaking checks
% %         a=find(H(i,:)>Hcr(i,:));    
% %         H(i,a)=Hcr(i,a);
% %         E(a)=Hcr(i,a).^2/16;
%         
% 
% 
% end;%end of wave propagation





        %[E]=wavepropagationstep(i,c,E,d,kwave,sigma,T,dx,stretchdx,Cbed,alpha,angle0,N,M,A,AW,p,ii,jj,periodic,Holateral);    
%         facspread=0.2;
%         Eo=E*(1-2*facspread);
%         Ei=E*facspread;
%         Ej=E*facspread;
%         
%         [Eo]=wavepropagationstep(i,c,cg,Eo,d,kwave,sigma,T,dx,stretchdx,Cbr,Cbed,alpha,angle0,N,M,A,AW,p,ii,jj,periodic,Holateral*sqrt(1-2*facspread));  
%         [Ei]=wavepropagationstep(i,c,cgi,Ei,di,kwavei,sigma,T,dx,stretchdx,Cbr,Cbed,alpha,angle0+spread,N,M,A,AWi,p,ii,jj,periodic,Holateral*sqrt(facspread)); 
%         [Ej]=wavepropagationstep(i,c,cgj,Ej,dj,kwavej,sigma,T,dx,stretchdx,Cbr,Cbed,alpha,angle0-spread,N,M,A,AWj,p,ii,jj,periodic,Holateral*sqrt(facspread)); 
%         
%         E=Eo+Ei+Ej;


       % [E]=wavepropagationstep(i,c,cgi,E,di,kwavei,sigma,T,dx,stretchdx,Cbr,Cbed,alpha,angle0+spread,N,M,A,AWi,p,ii,jj,periodic,Holateral); 
   

       
       
       
       
       
       
       
       
       
       
       
       
%        %Refraction
% if nrefrac>=1
%     a=find(dC(i,:)>0); %move from left to right (from i to o and from o to j)
%     if nrefrac>=4;dEi4i3=Ei4(a).*(1-1./(1+dCi4(i,a)/(spdi4/180*pi)));else;dEi4i3=0;end
%     if nrefrac>=3;dEi3i2=Ei3(a).*(1-1./(1+dCi3(i,a)/(spdi3/180*pi)));else;dEi3i2=0;end
%     if nrefrac>=2;dEi2i=Ei2(a).*(1-1./(1+dCi2(i,a)/(spdi2/180*pi)));else;dEi2i=0;end
%     dEio=Ei(a).*(1-1./(1+dCi(i,a)/(spdi/180*pi)));
%     dEoj=Eo(a).*(1-1./(1+dCj(i,a)/(spdj/180*pi)));
%     if nrefrac>=2;dEjj2=Ej(a).*(1-1./(1+dCj2(i,a)/(spdj2/180*pi)));else;dEjj2=0;end
%     if nrefrac>=3;dEj2j3=Ej2(a).*(1-1./(1+dCj3(i,a)/(spdj3/180*pi)));else;dEj2j3=0;end
%     if nrefrac>=4;dEj3j4=Ej3(a).*(1-1./(1+dCj4(i,a)/(spdj4/180*pi)));else;dEj3j4=0;end
% 
%     if nrefrac>=4;Ei4(a)=Ei4(a)-dEi4i3;end
%     if nrefrac>=3;Ei3(a)=Ei3(a)+dEi4i3-dEi3i2;end
%     if nrefrac>=2;Ei2(a)=Ei2(a)+dEi3i2-dEi2i;end
%     Ei(a)=Ei(a)+dEi2i-dEio;
%     Eo(a)=Eo(a)+dEio-dEoj;
%     Ej(a)=Ej(a)+dEoj-dEjj2;
%     if nrefrac>=2;Ej2(a)=Ej2(a)+dEjj2-dEj2j3;end
%     if nrefrac>=3;Ej3(a)=Ej3(a)+dEj2j3-dEj3j4;end
%     if nrefrac>=4;Ej4(a)=Ej4(a)+dEj3j4;end
% 
%     
%     a=find(dC(i,:)<0); %move from rigth to left
%     if nrefrac>=4;dEj4j3=Ej4(a).*(1-1./(1-dCj4(i,a)/(spdj4/180*pi)));else;dEj4j3=0;end
%     if nrefrac>=3;dEj3j2=Ej3(a).*(1-1./(1-dCj3(i,a)/(spdj3/180*pi)));else;dEj3j2=0;end
%     if nrefrac>=2;dEj2j=Ej2(a).*(1-1./(1-dCj2(i,a)/(spdj2/180*pi)));else;dEj2j=0;end
%     dEjo=Ej(a).*(1-1./(1-dCj(i,a)/(spdj/180*pi)));
%     dEoi=Eo(a).*(1-1./(1-dCi(i,a)/(spdi/180*pi)));
%     if nrefrac>=2;dEii2=Ei(a).*(1-1./(1-dCi2(i,a)/(spdi2/180*pi)));else;dEii2=0;end
%     if nrefrac>=3;dEi2i3=Ei2(a).*(1-1./(1-dCi3(i,a)/(spdi3/180*pi)));else;dEi2i3=0;end
%     if nrefrac>=4;dEi3i4=Ei3(a).*(1-1./(1-dCi4(i,a)/(spdi4/180*pi)));else;dEi3i4=0;end
% 
%     if nrefrac>=4;Ej4(a)=Ej4(a)-dEj4j3;end
%     if nrefrac>=3;Ej3(a)=Ej3(a)+dEj4j3-dEj3j2;end
%     if nrefrac>=2;Ej2(a)=Ej2(a)+dEj3j2-dEj2j;end
%     Ej(a)=Ej(a)+dEj2j-dEjo;
%     Eo(a)=Eo(a)+dEjo-dEoi;
%     Ei(a)=Ei(a)+dEoi-dEii2;
%     if nrefrac>=2;Ei2(a)=Ei2(a)+dEii2-dEi2i3;end  
%     if nrefrac>=3;Ei3(a)=Ei3(a)+dEi2i3-dEi3i4;end  
%     if nrefrac>=4;Ei4(a)=Ei4(a)+dEi3i4;end  
% end