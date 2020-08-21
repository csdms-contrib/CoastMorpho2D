function     [E,ALLE,dPW]=waveEVOLUTIONstep(EEE,E,ALLE,i,nrefrac,alpha,Cbr,Cbed,wavefrictionCollins,gridDIR,p,q,ii,jj,N,M,A,AW,AWup,d,dup,kwave,c,cg,cgup,bedfr,sigma,T,dx,periodic,Holateral,sigmatot,...
        ALLangle,ALLfacsprd,ALLspd,dCx,dCy);%  
    
    
        Eo=ALLE(1,:);if nrefrac>=1;Ei=ALLE(2,:);Ej=ALLE(3,:);end;if nrefrac>=2;Ei2=ALLE(4,:);Ej2=ALLE(5,:);end;
        if nrefrac>=3;Ei3=ALLE(6,:);Ej3=ALLE(7,:);end;if nrefrac>=4;Ei4=ALLE(8,:);Ej4=ALLE(9,:);end
        
        angle=ALLangle(1); 
        if nrefrac>=1
        anglei=ALLangle(2); 
        anglej=ALLangle(3);
        end
        if nrefrac>=2 
        anglei2=ALLangle(4); 
        anglej2=ALLangle(5); 
        end
        if nrefrac>=3
        anglei3=ALLangle(6); 
        anglej3=ALLangle(7); 
        end
        if nrefrac>=4
        anglei4=ALLangle(8); 
        anglej4=ALLangle(9); 
        end
        
        facsprdo=ALLfacsprd(1);
        if nrefrac>=1
        facsprdi=ALLfacsprd(2);
        facsprdj=ALLfacsprd(3);
        end
        if nrefrac>=2 
        facsprdi2=ALLfacsprd(4);
        facsprdj2=ALLfacsprd(5);
        end
        if nrefrac>=3 
        facsprdi3=ALLfacsprd(6);
        facsprdj3=ALLfacsprd(7);
        end
        if nrefrac>=4 
        facsprdi4=ALLfacsprd(8);
        facsprdj4=ALLfacsprd(9);
        end
        
        if nrefrac>=1
        spdi=ALLspd(1);
        spdj=ALLspd(2);
        end
        if nrefrac>=2 
        spdi2=ALLspd(3);
        spdj2=ALLspd(4);
        end
        if nrefrac>=3 
        spdi3=ALLspd(5);
        spdj3=ALLspd(6);
        end
        if nrefrac>=4 
        spdi4=ALLspd(7);
        spdj4=ALLspd(8);
        end
        

        
      
        
%Refraction
%Move from left to right (from i to o and from o to j)
%Move from rigth to left (from j to o and from o to i)
%INVERTED ORDER ON COTOBER 1st 2018

%     if nrefrac>=1;
%        tempEi=Ei;tempEj=Ej;tempEo=Eo;
%        dC=sin((angle+anglei)/2/180*pi)*dCx-cos((angle+anglei)/2/180*pi)*dCy;
%        a=find(dC>0);dEio=tempEi(a).*(1-1./(1+dC(a)/(spdi/180*pi)));Ei(a)=Ei(a)-dEio;Eo(a)=Eo(a)+dEio;
%        a=find(dC<0);dEoi=tempEo(a).*(1-1./(1-dC(a)/(spdi/180*pi)));Eo(a)=Eo(a)-dEoi;Ei(a)=Ei(a)+dEoi;
%        dC=sin((angle+anglej)/2/180*pi)*dCx-cos((angle+anglej)/2/180*pi)*dCy;
%        a=find(dC>0);dEoj=tempEo(a).*(1-1./(1+dC(a)/(spdj/180*pi)));Eo(a)=Eo(a)-dEoj;Ej(a)=Ej(a)+dEoj;
%        a=find(dC<0);dEjo=tempEj(a).*(1-1./(1-dC(a)/(spdj/180*pi)));Ej(a)=Ej(a)-dEjo;Eo(a)=Eo(a)+dEjo;
%     end
%     if nrefrac>=2;
%       tempEi=Ei;tempEj=Ej;tempEi2=Ei2;tempEj2=Ej2;
%       dC=sin((anglei+anglei2)/2/180*pi)*dCx-cos((anglei+anglei2)/2/180*pi)*dCy;
%       a=find(dC>0);dEi2i=tempEi2(a).*(1-1./(1+dC(a)/(spdi2/180*pi)));Ei2(a)=Ei2(a)-dEi2i;Ei(a)=Ei(a)+dEi2i;
%       a=find(dC<0);dEii2=tempEi(a).*(1-1./(1-dC(a)/(spdi2/180*pi)));Ei(a)=Ei(a)-dEii2;Ei2(a)=Ei2(a)+dEii2;
%       dC=sin((anglej+anglej2)/2/180*pi)*dCx-cos((anglej+anglej2)/2/180*pi)*dCy;
%       a=find(dC>0);dEjj2=tempEj(a).*(1-1./(1+dC(a)/(spdj2/180*pi)));Ej(a)=Ej(a)-dEjj2;Ej2(a)=Ej2(a)+dEjj2;
%       a=find(dC<0);dEj2j=tempEj2(a).*(1-1./(1-dC(a)/(spdj2/180*pi)));Ej2(a)=Ej2(a)-dEj2j;Ej(a)=Ej(a)+dEj2j;
%     end
%     if nrefrac>=3;
%       tempEi2=Ei2;tempEj2=Ej2;tempEi3=Ei3;tempEj3=Ej3;
%       dC=sin((anglei2+anglei3)/2/180*pi)*dCx-cos((anglei2+anglei3)/2/180*pi)*dCy;
%       a=find(dC>0);dEi3i2=tempEi3(a).*(1-1./(1+dC(a)/(spdi3/180*pi)));Ei3(a)=Ei3(a)-dEi3i2;Ei2(a)=Ei2(a)+dEi3i2;
%       a=find(dC<0);dEi2i3=tempEi2(a).*(1-1./(1-dC(a)/(spdi3/180*pi)));Ei2(a)=Ei2(a)-dEi2i3;Ei3(a)=Ei3(a)+dEi2i3;
%       dC=sin((anglei2+anglei3)/2/180*pi)*dCx-cos((anglei2+anglei3)/2/180*pi)*dCy;
%       a=find(dC>0);dEj2j3=tempEj2(a).*(1-1./(1+dC(a)/(spdj3/180*pi)));Ej2(a)=Ej2(a)-dEj2j3;Ej3(a)=Ej3(a)+dEj2j3;
%       a=find(dC<0);dEj3j2=tempEj3(a).*(1-1./(1-dC(a)/(spdj3/180*pi)));Ej3(a)=Ej3(a)-dEj3j2;Ej2(a)=Ej2(a)+dEj3j2;
%     end
%     if nrefrac>=4;
%       tempEi3=Ei3;tempEj3=Ej3;tempEi4=Ei4;tempEj4=Ej4;
%       dC=sin((anglei3+anglei4)/2/180*pi)*dCx-cos((anglei3+anglei4)/2/180*pi)*dCy;
%       a=find(dC>0);dEi4i3=tempEi4(a).*(1-1./(1+dC(a)/(spdi4/180*pi)));Ei4(a)=Ei4(a)-dEi4i3;Ei3(a)=Ei3(a)+dEi4i3;
%       a=find(dC<0);dEi3i4=tempEi3(a).*(1-1./(1-dC(a)/(spdi4/180*pi)));Ei3(a)=Ei3(a)-dEi3i4;Ei4(a)=Ei4(a)+dEi3i4;
%       dC=sin((anglei3+anglei4)/2/180*pi)*dCx-cos((anglei3+anglei4)/2/180*pi)*dCy;
%       a=find(dC>0);dEj3j4=tempEj3(a).*(1-1./(1+dC(a)/(spdj4/180*pi)));Ej3(a)=Ej3(a)-dEj3j4;Ej4(a)=Ej4(a)+dEj3j4;
%       a=find(dC<0);dEj4j3=tempEj4(a).*(1-1./(1-dC(a)/(spdj4/180*pi)));Ej4(a)=Ej4(a)-dEj4j3;Ej3(a)=Ej3(a)+dEj4j3;
%     end
%         
        
     

    if nrefrac>=1;
       tempEi=Ei;tempEj=Ej;tempEo=Eo;
       dC=sin((angle+anglei)/2/180*pi)*dCx-cos((angle+anglei)/2/180*pi)*dCy;
       a=find(dC>0);dEio=tempEi(a).*(1-1./(1+dC(a)/(spdi/180*pi)));Ei(a)=Ei(a)-dEio;Eo(a)=Eo(a)+dEio;
       a=find(dC<0);dEoi=tempEo(a).*(1-1./(1-dC(a)/(spdi/180*pi)));Eo(a)=Eo(a)-dEoi;Ei(a)=Ei(a)+dEoi;
       dC=sin((angle+anglej)/2/180*pi)*dCx-cos((angle+anglej)/2/180*pi)*dCy;
       a=find(dC>0);dEoj=tempEo(a).*(1-1./(1+dC(a)/(spdj/180*pi)));Eo(a)=Eo(a)-dEoj;Ej(a)=Ej(a)+dEoj;
       a=find(dC<0);dEjo=tempEj(a).*(1-1./(1-dC(a)/(spdj/180*pi)));Ej(a)=Ej(a)-dEjo;Eo(a)=Eo(a)+dEjo;
    end
    if nrefrac>=2;
      tempEi=Ei;tempEj=Ej;tempEi2=Ei2;tempEj2=Ej2;
      dC=sin((anglei+anglei2)/2/180*pi)*dCx-cos((anglei+anglei2)/2/180*pi)*dCy;
      a=find(dC>0);dEi2i=tempEi2(a).*(1-1./(1+dC(a)/(spdi2/180*pi)));Ei2(a)=Ei2(a)-dEi2i;Ei(a)=Ei(a)+dEi2i;
      a=find(dC<0);dEii2=tempEi(a).*(1-1./(1-dC(a)/(spdi2/180*pi)));Ei(a)=Ei(a)-dEii2;Ei2(a)=Ei2(a)+dEii2;
      dC=sin((anglej+anglej2)/2/180*pi)*dCx-cos((anglej+anglej2)/2/180*pi)*dCy;
      a=find(dC>0);dEjj2=tempEj(a).*(1-1./(1+dC(a)/(spdj2/180*pi)));Ej(a)=Ej(a)-dEjj2;Ej2(a)=Ej2(a)+dEjj2;
      a=find(dC<0);dEj2j=tempEj2(a).*(1-1./(1-dC(a)/(spdj2/180*pi)));Ej2(a)=Ej2(a)-dEj2j;Ej(a)=Ej(a)+dEj2j;
    end
    if nrefrac>=3;
      tempEi2=Ei2;tempEj2=Ej2;tempEi3=Ei3;tempEj3=Ej3;
      dC=sin((anglei2+anglei3)/2/180*pi)*dCx-cos((anglei2+anglei3)/2/180*pi)*dCy;
      a=find(dC>0);dEi3i2=tempEi3(a).*(1-1./(1+dC(a)/(spdi3/180*pi)));Ei3(a)=Ei3(a)-dEi3i2;Ei2(a)=Ei2(a)+dEi3i2;
      a=find(dC<0);dEi2i3=tempEi2(a).*(1-1./(1-dC(a)/(spdi3/180*pi)));Ei2(a)=Ei2(a)-dEi2i3;Ei3(a)=Ei3(a)+dEi2i3;
      dC=sin((anglei2+anglei3)/2/180*pi)*dCx-cos((anglei2+anglei3)/2/180*pi)*dCy;
      a=find(dC>0);dEj2j3=tempEj2(a).*(1-1./(1+dC(a)/(spdj3/180*pi)));Ej2(a)=Ej2(a)-dEj2j3;Ej3(a)=Ej3(a)+dEj2j3;
      a=find(dC<0);dEj3j2=tempEj3(a).*(1-1./(1-dC(a)/(spdj3/180*pi)));Ej3(a)=Ej3(a)-dEj3j2;Ej2(a)=Ej2(a)+dEj3j2;
    end
    if nrefrac>=4;
      tempEi3=Ei3;tempEj3=Ej3;tempEi4=Ei4;tempEj4=Ej4;
      dC=sin((anglei3+anglei4)/2/180*pi)*dCx-cos((anglei3+anglei4)/2/180*pi)*dCy;
      a=find(dC>0);dEi4i3=tempEi4(a).*(1-1./(1+dC(a)/(spdi4/180*pi)));Ei4(a)=Ei4(a)-dEi4i3;Ei3(a)=Ei3(a)+dEi4i3;
      a=find(dC<0);dEi3i4=tempEi3(a).*(1-1./(1-dC(a)/(spdi4/180*pi)));Ei3(a)=Ei3(a)-dEi3i4;Ei4(a)=Ei4(a)+dEi3i4;
      dC=sin((anglei3+anglei4)/2/180*pi)*dCx-cos((anglei3+anglei4)/2/180*pi)*dCy;
      a=find(dC>0);dEj3j4=tempEj3(a).*(1-1./(1+dC(a)/(spdj4/180*pi)));Ej3(a)=Ej3(a)-dEj3j4;Ej4(a)=Ej4(a)+dEj3j4;
      a=find(dC<0);dEj4j3=tempEj4(a).*(1-1./(1-dC(a)/(spdj4/180*pi)));Ej4(a)=Ej4(a)-dEj4j3;Ej3(a)=Ej3(a)+dEj4j3;
    end
        
        


    
    %%%SEGE  MRD  FATTOEE @ FRICTION OR NOT????
    
    %THIS IS THE GOOD ONE< FOR DELAWAREWE
    %bedfr=bedfr./Cbed.*(0.015*9.81.*(pi*(4*sqrt(E)/sqrt(2))./(T.*sinh(kwave.*max(d,0.1))))); %the sqrt(2) is to get Urms, not Upeak
    if wavefrictionCollins==1   
    bedfr=bedfr.*(1*0.015*9.81.*(pi*(4*sqrt(E)/sqrt(2))./(T.*sinh(kwave.*max(d,0.1))))); %the sqrt(2) is to get Urms, not Upeak
    else
    bedfr=bedfr*Cbed;    
    end
    
    %The factor 2 in front of 0.015.. Smbra che serva per Delaware, ma non
    %per VCR...??? (March 11th 2019)
    
    
    
    
    %bedfrCORR=1./Cbed.*(1*0.015*9.81.*max(0,(pi*(4*sqrt(E/Ejonswap)/sqrt(2))./T./sinh(kwave.*d)))); %the sqrt(2) is to get Urms, not Upeak
    %bedfr=bedfr./Cbed.*(1*0.015*9.81.*max(0,(pi*(4*sqrt(E)/sqrt(2))./min(10000,T.*sinh(kwave)))))*0; %the sqrt(2) is to get Urms, not Upeak
    %for the first two, cite Cavaleri, L., Malanotte-Rizzoli, P., 1981. Wind wave predictionin shallow water: theory and application. J
    %Energy balance of wind waves as a function of thebottom friction formulationR. Padilla-Hernandez) ´ , J. Monbaliu
    
    
    %anglei=angle;anglej=angle;
    %anglei2=angle;anglej2=angle;
    %anglei3=angle;anglej3=angle;
    %anglei4=angle;anglej4=angle;
%Propagate each energy according to its direction
    [Eo,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Eo,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,angle,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdo),sigmatot); 
    dPW=dPWs;
    if nrefrac>=1;
    [Ei,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ei,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglei,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdi),sigmatot);     
    dPW=dPW+dPWs;
    [Ej,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ej,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglej,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdj),sigmatot); 
    dPW=dPW+dPWs;
    end
    
    if nrefrac>=2;
    [Ei2,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ei2,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglei2,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdi2),sigmatot); 
    dPW=dPW+dPWs;
    [Ej2,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ej2,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglej2,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdj2),sigmatot); 
    dPW=dPW+dPWs;
    end
    
    if nrefrac>=3;
    [Ei3,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ei3,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglei3,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdi3),sigmatot); 
    dPW=dPW+dPWs;
    [Ej3,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ej3,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglej3,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdj3),sigmatot); 
    dPW=dPW+dPWs;
    end
    
    if nrefrac>=4;
    [Ei4,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ei4,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglei4,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdi4),sigmatot); 
    dPW=dPW+dPWs;
    [Ej4,dPWs]=wavepropagationONLYstepTEMPi(gridDIR,Ej4,E,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,anglej4,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral*sqrt(facsprdj4),sigmatot); 
    dPW=dPW+dPWs;
    end
    
    %Sum all the energies
%     E=Eo; 
%     if nrefrac>=1;    E=E+Ei+Ej;         end
%     if nrefrac>=2;    E=E+Ei2+Ej2;         end
%     if nrefrac>=3;    E=E+Ei3+Ej3;        end
%     if nrefrac>=4;    E=E+Ei4+Ej4;        end
     
    if nrefrac==0;    E=Eo;                   end
    if nrefrac==1;    E=Eo+Ei+Ej;            end
    if nrefrac==2;    E=Eo+Ei+Ej+Ei2+Ej2;          end
    if nrefrac==3;    E=Eo+Ei+Ej+Ei2+Ej2+Ei3+Ej3;          end
    if nrefrac==4;    E=Eo+Ei+Ej+Ei2+Ej2+Ei3+Ej3+Ei4+Ej4;          end
       

   
    
    
if nrefrac==0
ALLE=Eo;
elseif nrefrac==1
ALLE=[Eo;Ei;Ej];
elseif nrefrac==2
ALLE=[Eo;Ei;Ej;Ei2;Ej2];
elseif nrefrac==3
ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3];
elseif nrefrac==4
ALLE=[Eo;Ei;Ej;Ei2;Ej2;Ei3;Ej3;Ei4;Ej4];
end