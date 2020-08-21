function [z,FQsW_L,FQsW_R,longshore,angleshore]=bedevolutionDIV(deltaPW,U,fTide,A,AW,z,lev,wl,ho,Yreduction,N,M,dt,dx,Trange,Qs,rbulk,alpha,PWfactorlongshore,computewave,QsWslope,QsWon,angleswell,waveANGLE,Active,periodic,optionBC,gridDIR,FQsW_L,FQsW_R);
%NOTES
%remember in the min(h(p,0.1). These are "explicit term. You need to use
%the old values when re-evaluatign the fluxes

% B=A*0;
% B(:,1)=1;%the incomign boundary is the LEFT boundary
% B(:,end)=-1;%the incomign boundary is the RIGHT boundary

%ADDED ON MAY 2019, reduces computation when you have a lot of upland
A(Active==0)=0; %eliminate the cells in which the water depth is too small

%to transform to volumetric discharge
Qs=Qs/rbulk;
QsWslope=QsWslope/rbulk;
QsWon=QsWon/rbulk;


%Longshore transprot by radiation stresses (deltaPW is the cross-shore %gradient in wave power)
longshoretransportcoefficient=0.39/9.81*2.65/(2.65-1.030)*(24*3600)/rbulk; 
Qsradstress=longshoretransportcoefficient*deltaPW;

%downslope component of the rad stress
QsWslope_longshore=Qsradstress*15;


h=-z+wl;%water depth for morphology
hold=h;
D=(alpha*3600*24*Qs)*dt/(dx^2);

%REDUCTION FOR SED THICKNESS
QsWslope=QsWslope.*Yreduction;
QsWon=QsWon.*Yreduction;
Qsradstress=Qsradstress.*Yreduction;


%compute the shorelien orientation at every location
      MASK=h>0.2;
      [FX,FY] = gradient(-z);
      angleshore=atan2(FX,FY)*180/pi;
      angleshore(angleshore>80)=80;angleshore(angleshore<-80)=-80;

waveANGLElongshore=angleswell-angleshore;
waveANGLElongshore(waveANGLElongshore>80)=80;waveANGLElongshore(waveANGLElongshore<-80)=-80;

Qsradstress(waveANGLElongshore>=80 | waveANGLElongshore<=-80)=0;
QsWslope_longshore(waveANGLElongshore>=80 | waveANGLElongshore<=-80)=0;



%DO NOT MAKE DOWNSLOPE WITHIN BOUNDARY CELLS OR TO/FROM BOUNDARY CELLS
p=find(A==1);%exclude the NOLAND CELLS (A==0)| A==2
G=0*A;NN=length(p);G(p)=[1:NN];
rhs=h(p);
i=[];j=[];s=[];S=0*G;
for k = [N -N 1 -1]  % is left to right along the profile. If you take off N-N -> no coupling along-shore
    
    %boundary cells
    if periodic==0
    [a,q]=excludeboundarycell(k,N,M,p);%NOGRADIENT
    elseif periodic==1;
    [a,q]=periodicY(k,N,M,p); %for the long-shore
    end
    a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2


    valueC=(D(p(a))+D(q(a)))/2;
    
    
    valueWS=0*p(a); 
    valueWA=0*p(a);    
    Fout=0*A(p(a))+1;
    Fin=0*A(p(a))+1;
    if computewave==1;
        

        QsradstressI=(Qsradstress(p(a))+Qsradstress(q(a)))/2./max(1,h(p(a)))*dt/dx;
     
        valueWS=valueWS+max(QsWslope(p(a)),QsWslope(q(a)))*dt/dx^2;
       
        Qalongwave=0.2*QsWon(q(a))./max(0.1,h(p(a)))*dt/dx;%the 0.2 is the calibration coefficient after comparison with Dean'profle
        
     
        pry=cos(waveANGLE(p(a))/180*pi);
        prx=sin(abs(waveANGLE(p(a)))/180*pi);
                
        prylongshore=cos(waveANGLElongshore(p(a))/180*pi);%
        prxlongshore=sin(abs(waveANGLElongshore(p(a)))/180*pi);
        
        
        %Swell DOWNSLOPE flux due to BREAKING
        valueWS=valueWS+max(QsWslope_longshore(p(a)),QsWslope_longshore(q(a)))*dt/dx^2; 
         
   %CROSS SHORE by projection of alongwave
            if k==gridDIR;  %find the updrift direction!!!)
                valueWA=valueWA +pry.*Qalongwave;end %use the cell to the right!! UPWIND        
  
   %ALONG SHORE by projection of alongwave
            if (k==-N);aaa=find(waveANGLE(p(a))<0 );valueWA(aaa)=valueWA(aaa) +prx(aaa).*Qalongwave(aaa) ;end %use the cell to the right!! UPWIND
            if (k==N); aaa=find(waveANGLE(p(a))>=0);valueWA(aaa)=valueWA(aaa) +prx(aaa).*Qalongwave(aaa) ;end %use the cell to the right!! UPWIND            
                     
   %LONGSHORE by RADIATION STRESS
            if (k==-N);aaa=find(waveANGLElongshore(p(a))<0 );valueWA(aaa)=valueWA(aaa) +cos(angleshore(p(a(aaa)))/180*pi).*prylongshore(aaa).*prxlongshore(aaa).*QsradstressI(aaa);end %use the cell to the right!! UPWIND
            if (k==N); aaa=find(waveANGLElongshore(p(a))>=0);valueWA(aaa)=valueWA(aaa) +cos(angleshore(p(a(aaa)))/180*pi).*prylongshore(aaa).*prxlongshore(aaa).*QsradstressI(aaa);end %use the cell to the right!! UPWIND            
        %the first cos is to project the transport in the y direction
        %then there is a cos and sin as for the CERC formula
 
          if periodic==0
             if (k==-N);aaa=find(waveANGLE(p(a))<0);%length(aaa) %left to right  
             Fout(AW(q(a(aaa)))==1)=0;% sss=find(AW(q(a(aaa)))==1);length(sss)  %do not remove sediment from this cell
             Fin(AW(p(a(aaa)))==-1)=0;%  sss=find(AW(p(a(aaa)))==-1);length(sss)
             elseif (k==N);aaa=find(waveANGLE(p(a))>=0);%length(aaa) %right to left
             Fin(AW(p(a(aaa)))==1)=0; %sss=find(AW(p(a(aaa)))==1);length(sss)
             Fout(AW(q(a(aaa)))==-1)=0;%sss=find(AW(q(a(aaa)))==-1);length(sss)
             end
          end

    end
    
    
    %what enter from the cell (at the boundary already enter zero)
    S(p(a))=S(p(a))+(valueC+valueWS+valueWA.*Fin);
    %what exit from the cell. if double it scour too much. Need to put to zero to avoid scour at the inlet
    i=[i; G(q(a))];j=[j; G(p(a))]; s=[s; -(valueC+valueWS+valueWA.*Fout)]; 
end
%summary of the material that exits the cell
i=[i; G(p)]; j=[j; G(p)];s=[s; 1+S(p)];

%solve the matrix inversion
ds2=sparse(i,j,s);
P=ds2\rhs;


h(G>0)=full(P(G(G>0)));%rescale the matrix

%back to the bed elevation
z=wl-h;





























%%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!

if computewave==1;
if periodic==0
prx=sin(abs(waveANGLE)/180*pi);%*facLONGSHORE;
prylongshore=cos(waveANGLElongshore/180*pi);       
prxlongshore=sin(abs(waveANGLElongshore)/180*pi);        
   
%lateral boundary
for k = [N -N]
   
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2

%only get the cell that goes into a lateral cell
if (k==-N);    
    aaa=find(waveANGLE(p(a))<0);
    aL=a(AW(q(a(aaa)))==1);%length(aL)%find the left boundary
    FQsW_L=FQsW_L+sum(prx(p(aL)).*h(p(aL)).*0.2.*QsWon(q(aL))./max(0.1,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
    
    aR=a(AW(p(a(aaa)))==-1);%length(aR)%find the right boundary
    FQsW_R=FQsW_R-sum(prx(p(aR)).*h(p(aR)).*0.2.*QsWon(q(aR))./max(0.1,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));

 
end

if (k==N);    
    aaa=find(waveANGLE(p(a))>=0);
    aL=a(AW(p(a(aaa)))==1);%length(aL) %find the left boundary
    FQsW_L=FQsW_L-sum(prx(p(aL)).*h(p(aL)).*0.2.*QsWon(q(aL))./max(0.1,hold(p(aL)))*dt/dx);%%%%%%%%%%.*(hold(q(aL))>hold(p(aL))));
    
    aR=a(AW(q(a(aaa)))==-1);%length(aR)%find the right boundary
    FQsW_R=FQsW_R+sum(prx(p(aR)).*h(p(aR)).*0.2.*QsWon(q(aR))./max(0.1,hold(p(aR)))*dt/dx);%%%%%%%%%%.*(hold(q(aR))>hold(p(aR))));

 end


end
end
end
%%%%%%%%%%%%%%%%%%%








%This is not to balnce mass fluxes, it is just to make a plot of longshor
%transport for different pionts ont he laterl aboundary!!!
longshore=h*0;
if computewave==1;
prx=sin(abs(waveANGLE)/180*pi);%*facLONGSHORE;
%lateral boundary
for k = [N -N]
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2

%only get the cell that goes into a lateral cell
if (k==-N);    aaa=find(waveANGLE(p(a))<0);
    %longshore(p(aaa))=-prx(p(aaa)).*h(p(aaa)).*(QsWon(p(aaa))+QsWon(q(aaa)))/2./max(0.1,hold(p(aaa)));
end

if (k==N);    aaa=find(waveANGLE(p(a))>=0);
    %longshore(p(aaa))=prx(p(aaa)).*h(p(aaa)).*(QsWon(p(aaa))+QsWon(q(aaa)))/2./max(0.1,hold(p(aaa)));
end

end
end
%%%%%%%%%%%%%%%%%%%


% if angleswell>0;k=N;else;k=-N;end
% if periodic==0
% [a,q]=excludeboundarycell(k,N,M,p);
% elseif periodic==1;
% [a,q]=periodicY(k,N,M,p); %for the long-shore
% end
% %longshore=(prx(p).*h(p).*min(QsWon(p),QsWon(q))./max(0.1,hold(p))*dt/dx);
%         %prx=sin(abs(waveANGLE(p(a)))/180*pi)*facLONGSHORE;%./cos(waveANGLE(p(a))/180*pi);
% prx=sin(abs(waveANGLE)/180*pi)*facLONGSHORE;%./cos(waveANGLE(p(a))/180*pi);
% %longshore=prx.*h.*min(QsWon,QsWon);%./max(0.1,hold)*dt/dx);
% %QsWon(p)
% %QsWon(q)
% longshore(p)=prx(p).*h(p).*(QsWon(p)+QsWon(q))/2./max(0.1,hold(p));
% %longshore=prx.*QsWon;


%figure
%imagesc(longshore)
%pause


% %only get the cell that goes into a lateral cell
% if (angleswell<0 & k==-N);
%     aL=a(AW(q(a))==1);%find the left boundary
%     FQsW_L=FQsW_L+prx*sum(h(p(aL)).*min(QsWon(p(aL)),QsWon(q(aL)))./max(0.1,hold(p(aL)))*dt/dx);
%     
%     aR=a(AW(p(a))==-1);%find the right boundary
%     FQsW_R=FQsW_R-prx*sum(h(p(aR)).*min(QsWon(p(aR)),QsWon(q(aR)))./max(0.1,hold(p(aR)))*dt/dx);
% end
% 
% if (angleswell>0 & k==N);
%     aL=a(AW(p(a))==1); %find the left boundary
%     FQsW_L=FQsW_L-sum(prx(AW(p(a))==1).*h(p(aL)).*min(QsWon(p(aL)),QsWon(q(aL)))./max(0.1,hold(p(aL)))*dt/dx);
%     
%     aR=a(AW(q(a))==-1);%find the right boundary
%     FQsW_R=FQsW_R+sum(prx(AW(p(a))==1).*h(p(aR)).*min(QsWon(p(aR)),QsWon(q(aR)))./max(0.1,hold(p(aR)))*dt/dx);
%     %aL=a(AW(p(a))==1); %find the left boundary
%     %FQsW_L=FQsW_L-sum(prx.*h(p(aL)).*min(QsWon(p(aL)),QsWon(q(aL)))./max(0.1,hold(p(aL)))*dt/dx);
%     
%     %aR=a(AW(q(a))==-1);%find the right boundary
%     %FQsW_R=FQsW_R+sum(prx.*h(p(aR)).*min(QsWon(p(aR)),QsWon(q(aR)))./max(0.1,hold(p(aR)))*dt/dx);
% end