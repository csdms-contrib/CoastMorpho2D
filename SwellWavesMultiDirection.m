function [H,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,d,ho,hwSwell_lim,Tpeak,dx,periodic,angle,gridDIR,nrefrac,wavediffraction);
%OCIO YOU IMPOSED PERIODIC==1 in chsoign the q+k later on (foe Delaware %trick)

d(d<0 | A==0)=0.1;
%Energy balance of wind waves as a function of thebottom friction formulationR. Padilla-Hernandez) ´ , J. Monbaliu

%Parameters
alpha=50; %adimensional


spread=15;
spread2=30;
spread3=45;
spread4=60;


Holateral=Ho;  %this imposes the lateral boundary condition when the flag AW=2 is imposed. you can generaly choose either Ho or 0

    %DESCRIPTION
    %o is middle ray, 
    %i is left ray (angle turnerd clockwise, positive, increase the angle, whcih means will go toward lower left (more high angle).
    %j is right ray

% %condtions on the boundaries (directinoal spding imposed at the b.c.)   
%     facsprdo=1;
%     facsprdi=0.;
%     facsprdj=0.;
%     facsprdi2=0.;
%     facsprdj2=0.;
%     facsprdi3=0;
%     facsprdj3=0;
%     facsprdi4=0;
%     facsprdj4=0;
%%condtions on the boundaries (directinoal spding imposed at the b.c.)   
    facsprdo=0.6;
    facsprdi=0.1;
    facsprdj=0.1;
    facsprdi2=0.05;
    facsprdj2=0.05;
    facsprdi3=0.05;
    facsprdj3=0.05;
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



%kwave=0*d+1;kwave(d>0.2)=wavekgiulio(1/Tpeak,d(d>0.2));%kwave=wavek(1/Tp_swell,hwave);^Tpeak(j)
kwave=0*d+1;kwave(d>0.2)=wavek(1/Tpeak,d(d>0.2));%kwave=wavek(1/Tp_swell,hwave);^Tpeak(j)

Hcr=Cbr*d;
Hcr(d<=0.1)=0;
Hcr(d<=hwSwell_lim)=0;  

%IMPORTANT CHNAGE SEPT 2018!!!!!
Hcr(d<hwSwell_lim)=0;
    

sigma=2*pi/Tpeak;
c=(2*pi/Tpeak)./kwave;
cg=0.5*c.*(1+2*kwave.*d./sinh(2*kwave.*d));cg(d<=0.1 | cg<=0)=1; 
bedfr=1*(sigma./(9.81*sinh(kwave.*d))).^2;bedfr(d<=0.1)=1;
%NOTE: The Cbed or Collins will be multipled for later ---- bedfr=Cbed*(sigma./(9.81*sinh(kwave.*d))).^2;bedfr(d<=0.1)=1;


%Refraction coefficients
%if nrefrac>=1  %
%Note from giulio: there was a2 *2 in the past, a 0.5 around Dec 2018, and a
%1 in March 2019 (after taking off the 2 in the bed friction of Coolins)
dCx=([c(2:end,:); c(end,:)]-[c(1,:); c(1:end-1,:)])/2./c.*(d>1)*1; %NS
dCy=([c(:,2:end) c(:,end)]-[c(:,1) c(:,1:end-1)])/2./c.*(d>1)*1;   %EW
%end



H=zeros(N,M);H(end,:)=Ho;
waveANGLE=d*NaN;%just to initialize
wavePERIOD=d*0;%just to initialize
PWswell=d*0;   
deltaPW=d*0;   
    


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
E=H(end,:).^2/16;
    
%Partition the wave heigth in the boundary condtions    

    Eo=E*facsprdo;
    if nrefrac>=1;    Ei=E*facsprdi;      Ej=E*facsprdj;    end
    if nrefrac>=2;    Ei2=E*facsprdi2;    Ej2=E*facsprdj2;    end
    if nrefrac>=3;    Ei3=E*facsprdi3;    Ej3=E*facsprdj3;    end;
    if nrefrac>=4;    Ei4=E*facsprdi4;    Ej4=E*facsprdj4;    end
if nrefrac==0;ALLE=Eo;elseif nrefrac==1;ALLE=[Eo;Ei;Ej];elseif nrefrac==2;ALLE=[Eo;Ei;Ej;Ei2;Ej2];elseif nrefrac==3;ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4;ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACTUAL LOOP from the wave boundary and up following the waves  
for i=N-1:-1:1  
    
    
%This slighlty change the results on the boundaries... very small    
%sigmatot=2*pi*(max(0,E1./E)/Tpeak);sigmatot(E==0)=0;   
%wavePERIOD(i,sigmatot>0)=2*pi./sigmatot(sigmatot>0);
%in theory you should put sigmatotp == 0 where E is zero!!!
sigmatot=2*pi/Tpeak; 
wavePERIOD(i,:)=Tpeak;

%PWswell(i,:)=cg(i+1,:)*1030*9.8.*E;%maybe this ok ?
PWswell(i+1,:)=cg(i+1,:)*1030*9.8.*E;

%[E,ALLE]=waveEVOLUTIONstep(E,E,ALLE,i,nrefrac,alpha,Cbr,Cbed,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave(i,:),c,cg(i,:),cg(i+1,:),bedfr(i,:),sigma,Tpeak,dx,periodic,Holateral,sigmatot,...
%       ALLangle,ALLfacsprd,ALLspd,dCx(i,:),dCy(i,:));


%%THA GOOD
  [E,ALLE,dPW]=waveEVOLUTIONstep(E,E,ALLE,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave(i,:),c,cg(i,:),cg(i+1,:),bedfr(i,:),sigma,Tpeak,dx,periodic,Holateral,sigmatot,...
        ALLangle,ALLfacsprd,ALLspd,dCx(i+1,:),dCy(i+1,:));

 % [E,ALLE,dPW]=waveEVOLUTIONstep(E,E,ALLE,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW(i,:),AW(i+1,:),d(i,:),d(i+1,:),kwave(i+1,:),c,cg(i,:),cg(i+1,:),bedfr(i+1,:),sigma,Tpeak,dx,periodic,Holateral,sigmatot,...
%        ALLangle,ALLfacsprd,ALLspd,dCx(i+1,:),dCy(i+1,:));



%deltaPW(i,:)=dPW*1030*9.8;

ffE=max(0,E./E);    

     if nrefrac>=1;
        %%%asses ratio of energy from different directions
        f1Eo=max(0,ALLE(1,:)./E);     f1Ei=max(0,ALLE(2,:)./E);          f1Ej=max(0,ALLE(3,:)./E); 
        if nrefrac>=2;        f1Ei2=max(0,ALLE(4,:)./E);        f1Ej2=max(0,ALLE(5,:)./E);         end
        if nrefrac>=3;        f1Ei3=max(0,ALLE(6,:)./E);        f1Ej3=max(0,ALLE(7,:)./E);         end
        if nrefrac>=4;        f1Ei4=max(0,ALLE(8,:)./E);        f1Ej4=max(0,ALLE(9,:)./E);         end
     end
    
        wANGtot=0;    
        wANG=angle*ones(1,M);
        if nrefrac>=1;        wANG=angle*f1Eo+anglei*f1Ei+anglej*f1Ej;        end
        if nrefrac>=2;        wANG=wANG+anglei2*f1Ei2+anglej2*f1Ej2;        end
        if nrefrac>=3;        wANG=wANG+anglei3*f1Ei3+anglej3*f1Ej3;        end
        if nrefrac>=4;        wANG=wANG+anglei4*f1Ei4+anglej4*f1Ej4;        end
        wANGtot=wANGtot+wANG.*ffE;
    
        waveANGLE(i,:)=wANGtot;  
        
    %Diffract/Diffuse
    if wavediffraction==1;
    Mc=c(i,:);Mcg=cg(i,:);
    [E]=wavediffusionstep(E,i,Mc,Mcg,d(i,:),ho(i,:),sigmatot,dx,angle,Cbr,Cbed,alpha,angle,N,M,A,AW(i,:),AW(i+1,:),p,ii,jj,periodic,Holateral); 
    end
    
    
    
    %Recalculate the wave height and set breaking criterion check
    H(i,:)=4*sqrt(max(0,E));
    a=find(H(i,:)>Hcr(i,:));%waves will break!
    %deltaPW(i+1,a)=(H(i,a).^2-Hcr(i,a).^2)/dx*1030*9.8/16.*cg(i,a);
    H(i,a)=Hcr(i,H(i,:)>Hcr(i,:));
    E(a)=Hcr(i,a).^2/16; 
    
    %Breaking for longshore transport
    a=find(H(i+1,:)>Hcr(i,:));%waves will break!
    %deltaPW(i+1,a)=max(0,H(i+1,a).^2.*cg(i+1,a)-Hcr(i,a).^2.*cg(i,a))/dx*1030*9.8/16; 
    deltaPW(i+1,a)=max(0,H(i+1,a).^2.*cg(i+1,a)-Hcr(i,a).^2.*cg(i,a))/dx*1030*9.8/16; 
    %This is the radiation stress with the right units!!!
    
    %E=ffE.*E;
    ALLE(1,:)=E;
    if nrefrac>=1;        ALLE(1,:)=E.*f1Eo; ALLE(2,:)=E.*f1Ei;         ALLE(3,:)=E.*f1Ej;        end
    if nrefrac>=2;        ALLE(4,:)=E.*f1Ei2;        ALLE(5,:)=E.*f1Ej2;       end
    if nrefrac>=3;        ALLE(6,:)=E.*f1Ei3;        ALLE(7,:)=E.*f1Ej3;       end
    if nrefrac>=4;        ALLE(8,:)=E.*f1Ei4;        ALLE(9,:)=E.*f1Ej4;       end
  
end;%end of wave evolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uwave=(pi*H./Tpeak./sinh(kwave.*d));Uwave(d<=0.1)=0;

%flip the domain back
if gridDIR==-1; 
    H=H(end:-1:1,:);
    wavePERIOD=wavePERIOD(end:-1:1,:);
    waveANGLE=waveANGLE(end:-1:1,:);
    PWswell=PWswell(end:-1:1,:);
    Uwave=Uwave(end:-1:1,:);
    kwave=kwave(end:-1:1,:);
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