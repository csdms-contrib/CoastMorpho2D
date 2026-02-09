function [z,FQsW_L,FQsW_R,longshore,angleshore]=bedevolutionDIV(deltaPW,U,fTide,A,AW,z,Zlev,wl,ho,Y,Yreduction,VEG,N,M,dt,dx,Trange,crSAND,hlimbankro,hlimbankr,drybankbeta,Qs,reduceSANDbankVEG,rbulk,hwSwelltransport_lim,computewave,QsWslope,QsWon,angleswell,waveANGLE,Active,periodic,wavetransportzerolateralgradient,gridDIR,FQsW_L,FQsW_R);
%NOTES
%remember in the min(h(p,0.1). These are "explicit term. You need to use
%the old values when re-evaluatign the fluxes

% AW(:,1)=1;%the incomign boundary is the LEFT boundary
% AW(:,end)=-1;%the incomign boundary is the RIGHT boundary


%Reduces computation when you have a lot of upland
A(Active==0)=0; %eliminate the cells in which the water depth is too small


h=wl-z;%water depth for morphology
hold=h;

%figure;plot(wl);pause

QsWon=0.2*QsWon;%0.2 %This is really the alongwave, not the onshore (the onshore is Qswon projected on x)
QsWslope=1*QsWslope;

hminalongwave=0.1;%this is just numeric, shoudl be very small

%to transform to volumetric discharge
Qs=Qs/rbulk;
QsWslope=QsWslope/rbulk;
QsWon=QsWon/rbulk;

hminlongshore=2; %this can be larger, on the order of meters

%longshoretransportcoefficient=0.2/9.81*2.65/(2.65-1.030)*(24*3600)/rbulk; 
longshoretransportcoefficient=0.5/9.81*2.65/(2.65-1.030)*(24*3600)/rbulk; 

%deltaPW is the cross-shore gradient in wave power!!!
Qsradstress=deltaPW*longshoretransportcoefficient;

%Downslope component of the longshore transport caused by radiation stress wave breakoing)
QsWslope_longshore=deltaPW*longshoretransportcoefficient*15;%    *2;%*2;%2.5;%4;  %needed to avoid le buche!!!

%D=(3600*24*Qs)*dt/(dx^2);
D=(3600*24*Qs +crSAND)*dt/(dx^2);



%REDUCTION FOR SED THICKNESS
QsWon=QsWon.*Yreduction;
QsWslope=QsWslope.*Yreduction;
QsWslope_longshore=QsWslope_longshore.*Yreduction;
%%% Qs=Qs.*Yreduction;THIS IS DONE LATER IN THE CODE


%SHORELINE ORIENTATION
%MASK=h>0.2;%A~=0;%
%zz=diffusebedforshoreangle(z,A,MASK,0.01,dx,periodic);%zz=diffusefetch(A*0+1,z,10,dx); %F(F<=Fetchlim)=0;
%zz=z;zz(A==0)=NaN;
[FX,FY] = gradient(-z);
%FX=diffusebedforshoreangle(FX,A,MASK,1,dx,periodic);%zz=diffusefetch(A*0+1,z,10,dx); %F(F<=Fetchlim)=0;
%FY=diffusebedforshoreangle(FY,A,MASK,1,dx,periodic);%zz=diffusefetch(A*0+1,z,10,dx); %F(F<=Fetchlim)=0;
angleshore=atan2(FX,FY)*180/pi;%angleshore(A==0)=NaN;
%angleshoreI=angleshore;  

waveANGLElongshore=angleswell-angleshore;

Qsradstress=Qsradstress.*(angleshore<80 & angleshore>-80).*(waveANGLElongshore<80 & waveANGLElongshore>-80);
          
%figure;imagesc(z);pause

limiteforlatertranp=(h>hwSwelltransport_lim);
                      
pry=cos(waveANGLE/180*pi);
prx=sin(abs(waveANGLE)/180*pi).*limiteforlatertranp;
               
prylongshore=cos(waveANGLElongshore/180*pi);         
prxlongshore=sin(abs(waveANGLElongshore)/180*pi).*limiteforlatertranp;%.*(h(q(a))>h(p(a)));%.*slope;%(0.5+0.5*(h(q(a))>h(p(a))));%*facLONGSHORE  %ZIIIOOOOO        



%DO NOT MAKE CREEP WITHIN BOUNDARY CELLS OR TO/FROM BOUNDARY CELLS
p=find(A==1);%exclude the NOLAND CELLS (A==0)| A==2
G=0*A;NN=length(p);G(p)=[1:NN];
rhs=h(p);i=[];j=[];s=[];S=0*G;
for k = [N -N 1 -1]  % is left to right along the profile. If you take off N-N -> no coupling along-shore
    %boundary cells
    if periodic==0
    [a,q]=excludeboundarycell(k,N,M,p);%NOGRADIENT
    elseif periodic==1;
    [a,q]=periodicY(k,N,M,p); %for the long-shore
    end
    a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2
   
    %Current downslope flux Qs 
    %valueC=max(D(p(a)),D(q(a))).*(Yreduction(p(a))+Yreduction(q(a)))/2   .*(Y(p(a))>-0.5 & Y(q(a))>-0.5); 
 
    %USED FOR BARNSTA AND PIE PAPERS
    %valueC=max(D(p(a)),D(q(a))).*(Yreduction(p(a))+Yreduction(q(a)))/2   .*max( (z(p(a))>z(q(a))).*(Y(p(a))>-0.1) ,(z(q(a))>z(p(a))).*(Y(q(a))>-0.1));%  l      (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;
    
    %USED FOR NEW DELTA MODEL AFTER Jan 2024
    %acurally we dont continue using in the main code after March 2024
   
    
   %New
%    aUPWND=0.99;
%     valueC=(  (z(p(a))>z(q(a))).*(D(p(a))*aUPWND+D(q(a))*(1-aUPWND))   +(z(q(a))>z(p(a))).*(D(q(a))*aUPWND+D(p(a))*(1-aUPWND))   )...
%     .*max( (z(p(a))>z(q(a))).*Yreduction(p(a)) ,(z(q(a))>z(p(a))).*Yreduction(q(a)) );    %  (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;
% %     

%New 8 19 2024
%1)
%    ap=min(1,max(0,(h(p(a))-hDlo)/hDhi));ap=ap*(pDm-pDo)+pDo;
%    aq=min(1,max(0,(h(q(a))-hDlo)/hDhi));aq=aq*(pDm-pDo)+pDo;
%    valueC=(  (z(p(a))>z(q(a))).*(D(p(a)).*(1-ap) +D(q(a)).*ap)...
%            +(z(q(a))>z(p(a))).*(D(q(a)).*(1-aq) +D(p(a)).*aq)   )...
%   .*max( (z(p(a))>z(q(a))).*Yreduction(p(a)) ,(z(q(a))>z(p(a))).*Yreduction(q(a)) );% .*max(h(p(a)),h(q(a)))./((h(p(a))+h(q(a)))/2) ;    %  (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;
%    
% valueC=valueC.*(h(p(a))>0.1 & h(q(a))>0.1);

%2)
%valueC=max(D(p(a)),D(q(a)))...
%   .*max( (z(p(a))>z(q(a))).*Yreduction(p(a)) ,(z(q(a))>z(p(a))).*Yreduction(q(a)) );% .*max(h(p(a)),h(q(a)))./((h(p(a))+h(q(a)))/2) ;    %  (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;
   

%3)
%hlimbankro=0.01;%.1;%limited for dry bank erosion
%hlimbankr=0.5;%.1;%limited for dry bank erosion
%redu=1;% how much to reduce downslope when upslope cell is "dry", so swith mechanism of erosion.

%P.hlimbankro=0.01;%.1;%limited for dry bank erosion
%P.hlimbankr=0.5;%0.2;%.1;%limited for dry bank erosion
%P.drybankbeta=0.5;% how much to reduce downslope when upslope cell is "dry", so swith mechanism of erosion.%valueC=max( D(p(a)).*(1-redu*(h(q(a))<hlimbankr)) , D(q(a)).*(1-redu*(h(p(a))<hlimbankr)) )...

aq=(hlimbankr-h(q(a)))/(hlimbankr-hlimbankro);aq(h(q(a))>hlimbankr)=0;aq(h(q(a))<hlimbankro)=1;
ap=(hlimbankr-h(p(a)))/(hlimbankr-hlimbankro);ap(h(p(a))>hlimbankr)=0;ap(h(p(a))<hlimbankro)=1;
valueC=max( D(p(a)).*(1-(1-drybankbeta)*aq), D(q(a)).*(1-(1-drybankbeta)*ap))...
   .*max( (z(p(a))>z(q(a))).*Yreduction(p(a)) ,(z(q(a))>z(p(a))).*Yreduction(q(a)) );% .*max(h(p(a)),h(q(a)))./((h(p(a))+h(q(a)))/2) ;    %  (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;



   %BESTA USED THE MOST
   %valueC=(D(p(a))+D(q(a)))/2.*max( (z(p(a))>z(q(a))).*Yreduction(p(a)) ,(z(q(a))>z(p(a))).*Yreduction(q(a)) );%   .*max( (z(p(a))>z(q(a))).*(Y(p(a))>0) ,(z(q(a))>z(p(a))).*(Y(q(a))>0));%  l        (Y(p(a))>-0.5 & Y(q(a))>-0.5);   %valueC=(D(p(a))+D(q(a)))/2.*(Yreduction(p(a))+Yreduction(q(a)))/2;

   
   
   
   
   %erosion reduction at the vegeated bank
   valueC( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) )= ...
   reduceSANDbankVEG*valueC( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) );
%else
%    valueC=0*z(p(a));   
%end

% 
% %reduce downslope if dry cells
% hm=min(h(q(a)),h(p(a)));
% valueC(hm<0.2)=valueC(hm<0.2)*0.1;
%

    
    valueWS=0*p(a); 
    valueWA=0*p(a);  
    valueWL=0*p(a);     
    Fout=0*A(p(a))+1;Fin=0*A(p(a))+1;  
    FoutL=0*A(p(a))+1;FinL=0*A(p(a))+1;
    
    if computewave==1;   
        
    %Swell alongwave
    QsWonVALUE=QsWon(q(a))./max(hminalongwave,h(p(a)))*dt/dx;%*2;  %BETTER%<---- OCIO QUESTO *2 7pm april 24  2019
       

    %Radiation stress
    Qlonsh=(Qsradstress(p(a))+Qsradstress(q(a)))/2./max(hminlongshore,h(p(a)))*dt/dx;    %Qlonsh(isnan(angleshoreI))=0;      
      
   
    %Waves (Swell+Sea) DOWNSLOPE flux 
    valueWS=valueWS+(QsWslope(p(a))+QsWslope(q(a)))/2*dt/dx^2 .*(VEG(p(a))==0 & VEG(q(a))==0);    %.*min(Yreduction(p(a)),Yreduction(q(a)));% .*(VEG(p(a))==0 & VEG(q(a))==0);   

    %Swell DOWNSLOPE flux due to BREAKING
    valueWS=valueWS+(QsWslope_longshore(p(a))+QsWslope_longshore(q(a)))/2*dt/dx^2 ; %.*min(Yreduction(p(a)),Yreduction(q(a)));
         
    %CROSS SHORE
    if k==gridDIR;  %find the updrift direction
    valueWA=valueWA +pry(p(a)).*QsWonVALUE;
    end        
  
    %ALONGSHORE - USE UPWIND!
    if (k==-N);aaa=find(waveANGLE(p(a))<0 );valueWA(aaa)=valueWA(aaa) +prx(p(a(aaa))).*QsWonVALUE(aaa) ;end %use the cell to the right!! UPWIND
    if (k==N); aaa=find(waveANGLE(p(a))>=0);valueWA(aaa)=valueWA(aaa) +prx(p(a(aaa))).*QsWonVALUE(aaa) ;end %use the cell to the right!! UPWIND                     
            
    if (k==-N);aaa=find(waveANGLElongshore(p(a))<0 );valueWL(aaa)=valueWL(aaa) +cos(angleshore(p(a(aaa)))/180*pi).*prylongshore(p(a(aaa))).*prxlongshore(p(a(aaa))).*Qlonsh(aaa);end %use the cell to the right!! UPWIND
    if (k==N); aaa=find(waveANGLElongshore(p(a))>=0);valueWL(aaa)=valueWL(aaa) +cos(angleshore(p(a(aaa)))/180*pi).*prylongshore(p(a(aaa))).*prxlongshore(p(a(aaa))).*Qlonsh(aaa);end %use the cell to the right!! UPWIND            


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if periodic==0 & wavetransportzerolateralgradient==1
             if (k==-N)
             Fout(AW(q(a))==1 & waveANGLE(p(a))<0)=0;
             Fin(AW(p(a))==-1 & waveANGLE(p(a))<0)=0;
             elseif (k==N)
             Fout(AW(q(a))==-1 & waveANGLE(p(a))>=0)=0;
             Fin(AW(p(a))==1   & waveANGLE(p(a))>=0)=0;
             end  
             
             if (k==-N)
             FoutL(AW(q(a))==1 & waveANGLElongshore(p(a))<0)=0;
             FinL(AW(p(a))==-1 & waveANGLElongshore(p(a))<0)=0;
             elseif (k==N)
             FoutL(AW(q(a))==-1 & waveANGLElongshore(p(a))>=0)=0;
             FinL(AW(p(a))==1   & waveANGLElongshore(p(a))>=0)=0;
             end   
          end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end %End compute wave

    %what enter from the cell (at the boundary already enter zero)
    S(p(a))=S(p(a))+(valueC+valueWS+valueWA.*Fin+valueWL.*FinL);
    %what exit from the cell. if double it scour too much. Need to put to zero to avoid scour at the inlet
    i=[i; G(q(a))];j=[j; G(p(a))]; s=[s; -(valueC+valueWS+valueWA.*Fout+valueWL.*FoutL)]; 
end
%summary of the material that exits the cell
i=[i; G(p)]; j=[j; G(p)];s=[s; 1+S(p)];

ds2=sparse(i,j,s);P=ds2\rhs;%solve the matrix inversion

h(G>0)=full(P(G(G>0)));%rescale the matrix

%back to the bed elevation
z=wl-h;








                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!

if computewave==1;
if periodic==0 & wavetransportzerolateralgradient==1

%lateral boundary
for k = [N -N]
   
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2


%%%QSWon %only get the cell that goes into a lateral cell
if (k==-N);  
    aL=a(AW(q(a))==1 & waveANGLE(p(a))<0);%length(aL)%find the left boundary
    FQsW_L=FQsW_L+sum(prx(p(aL)).*h(p(aL)).*QsWon(q(aL))./max(hminalongwave,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
    
    aR=a(AW(p(a))==-1 & waveANGLE(p(a))<0);%length(aR)%find the right boundary
    FQsW_R=FQsW_R-sum(prx(p(aR)).*h(p(aR)).*QsWon(q(aR))./max(hminalongwave,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));
end
if (k==N);   
    aL=a(AW(p(a))==1  & waveANGLE(p(a))>=0);%length(aL) %find the left boundary
    FQsW_L=FQsW_L-sum(prx(p(aL)).*h(p(aL)).*QsWon(q(aL))./max(hminalongwave,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
    
    aR=a(AW(q(a))==-1 & waveANGLE(p(a))>=0);%length(aR)%find the right boundary
    FQsW_R=FQsW_R+sum(prx(p(aR)).*h(p(aR)).*QsWon(q(aR))./max(hminalongwave,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));
end


%%%LONGSHORE
if (k==-N) 
    aL=a(AW(q(a))==1 & waveANGLElongshore(p(a))<0);%length(aL)%find the left boundary
    FQsW_L=FQsW_L+sum(cos(angleshore(p((aL)))/180*pi).*prylongshore(p(aL)).*prxlongshore(p(aL))  .*h(p(aL)).*(Qsradstress(p(aL))+Qsradstress(q(aL)))/2./max(hminlongshore,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
    
    aR=a(AW(p(a))==-1 & waveANGLElongshore(p(a))<0);%length(aR)%find the right boundary
    FQsW_R=FQsW_R-sum(cos(angleshore(p((aR)))/180*pi).*prylongshore(p(aR)).*prxlongshore(p(aR))  .*h(p(aR)).*(Qsradstress(p(aR))+Qsradstress(q(aR)))/2./max(hminlongshore,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));
end
if (k==N)
    aL=a(AW(p(a))==1  & waveANGLElongshore(p(a))>=0);%length(aL) %find the left boundary
    FQsW_L=FQsW_L-sum(cos(angleshore(p((aL)))/180*pi).*prylongshore(p(aL)).*prxlongshore(p(aL))  .*h(p(aL)).*(Qsradstress(p(aL))+Qsradstress(q(aL)))/2./max(hminlongshore,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
   
    aR=a(AW(q(a))==-1 & waveANGLElongshore(p(a))>=0);%length(aR)%find the right boundary
    FQsW_R=FQsW_R+sum(cos(angleshore(p((aR)))/180*pi).*prylongshore(p(aR)).*prxlongshore(p(aR))  .*h(p(aR)).*(Qsradstress(p(aR))+Qsradstress(q(aR)))/2./max(hminlongshore,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));
end


end
end
end
%%%%%%%%%%%%%%%%%%%





%%%
longshore=NaN;
