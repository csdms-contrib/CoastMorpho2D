function [EmD,SSM,FLX,XXTide,XXRiver]=sedtran(numeric,d,A,SPCLcell,Do,DiffS,Zlev,h,ho,E,WS,dx,dt,rbulk,co,SeaSSCbelowLIMIT,ZlevcoLIMIT,Ux,Uy,fTide,Ttide,q1,qm1,qN,qmN,Cq1,CqN,periodic,computeriver,computetide,residualcurrents,kro,coRIVER,FLX,tracksedimentfluxes,XX);

Cqm1=[ Cq1(1,:); Cq1(1:end-1,:)];CqmN=[ CqN(:,1) CqN(:,1:end-1)];         


%%%%USED for ACE Basin, south carolina
%A(A==0)=2;

rivermouthfront=SPCLcell.rivermouthfront;

QmouthRiver=FLX(1);
QseaTide=FLX(2);
QseaRiver=FLX(3);
QmouthTide=FLX(4);
%consider the pond cells as A==0
%A(A==3)=1; %but do not update it!
%A(A==22)=2;

%A(ho<=0.001)=0; %ATTENTION: eliminate the cells in which the water depth is too small. They cannot evolve at all!!!! They will not erode nor accrete!!!!!!!!!!!!!!!!!!!!!!
A(rivermouthfront)=1;%the cells in front of the river mouth: let them erode if needed


%rivermouthfront
p = find(A>0);%exclude the NOLAND CELLS (A==0)
G=0*d;NN=length(p);G(p)=[1:NN];
rhs=E(p); %in the rhs there is already the additon of the erosion input
[N,M]=size(G);

%boundary conditions imposed SSC
if SeaSSCbelowLIMIT==1;
a=find(A==2);rhs(G(a))=co.*h(a).*(Zlev(a)<ZlevcoLIMIT);
else
a=find(A==2);rhs(G(a))=co.*h(a);
end

%iver and mud (if mud, you impose the SSC at the inlet)
if computeriver==1;
    for i=1:length(coRIVER)
    a=find(A==(10+i-1));rhs(G(a))=coRIVER(i).*h(a);%.*fTide(a);  
    end
end



Dxx=(Do*24*3600 +DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
Dyy=(Do*24*3600 +DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);

%Dxx(1:2,:)=0;Dyy(1:2,:)=0; %this one is the problem
%Dxx(end-1:end,:)=0;Dyy(end-1:end,:)=0;


% Dbase=(Do*24*3600)/(dx^2).*h;
% qx=(q1+qm1)/2;
% qy=(qN+qmN)/2;
% qq=sqrt(qx.^2+qy.^2);
% qx=abs(qx)./qq;
% qy=abs(qy)./qq;
% Dxx=Dbase.*qx;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
% Dyy=Dbase.*qy;%.*h;%.*(ho>kro);%.*(hgross>0.01);
% end


%Crazy test - SHOUDL NEVER BE USED!
% % % Dtot=sqrt(Dxx.^2+Dyy.^2);
% % % Dxx=Dtot;
% % % Dyy=Dtot;







%Dxx(A>=10 & A<=19)=0;
%Dyy(A>=10 & A<=19)=0; %no tidal flux at the river mouth
%Dxx(rivermouthfront)=0;Dyy(rivermouthfront)=0; %no tidal flux at the river mouth
%the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day

%Dxx=Dxx./fTide;
%Dyy=Dyy./fTide;

%teste remvoed jan 2025
%Dxx=Dxx.*max(0.01,fTide);
%Dyy=Dyy.*max(0.01,fTide);

%removed feb 2025
Dxx=Dxx.*fTide;
Dyy=Dyy.*fTide;





i=[];j=[];s=[];S=0*G;
%[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]  %go over the 4 directions for the gradients

%avoid to the the cells out of the domain (risk to make it periodic...)
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

a=a(A(q(a))==1 | A(q(a))==2 | (A(q(a))>=10 & A(q(a))<=19));%exlcude the translated cell that are NOLAND cells

%go over the extradiagonal direction for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;


DD=(D(p(a))+D(q(a)))/2;%.*min(fTide(p(a)),fTide(q(a)));%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
DD(isnan(DD))=0;%new Oct 2025
% if numeric==1
% DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
% elseif numeric==0
%DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
% end

Fin=0*DD;Fin(A(p(a))==1)=1; % (A(q(a))==1 | A(q(a))==10)=1; to conserve the mass a the river mouth = no input
Fout=0*DD;Fout(A(q(a))==1)=1; %needed not to affect the b.c. -> Do not choose 2 and 1p

%tidal dispersion component
value=DD./h(p(a));%./fTide(p(a));



%figure;imagesc(URy);pause
%river flow component
if computeriver==1
    %figure;imagesc(qN);pause

up=[];F=0;
if (k==N);UR=qN(p(a));up=find(UR>=0);F=UR(up);end %East-west
if (k==-N);UR=qmN(p(a));up=find(UR<0);F=-UR(up);end
if (k==1);UR=q1(p(a));up=find(UR>=0);F=UR(up);end  %North-south
if (k==-1);UR=qm1(p(a));up=find(UR<0);F=-UR(up);end
%value(up)=value(up)+F*3600*24/dx./h(p(a(up))).*min(fTide(p(a(up))),fTide(q(a(up)))).^1;%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;
%value(up)=value(up)+F*3600*24/dx./h(p(a(up))).*(fTide(p(a(up)))+fTide(q(a(up))))/2;%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;BEST AUG 5th
value(up)=value(up)+F*3600*24/dx./h(p(a(up))).*min(fTide(p(a(up))),fTide(q(a(up))));%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;
%value(up)=value(up)+F*3600*24/dx./h(p(a(up))).*min(fTide(p(a(up))),fTide(q(a(up))));
%.*min(h(p(a(up))),h(q(a(up))))./(0.5*(h(p(a(up)))+h(q(a(up)))))
% 
% up=[];F=0;
% if (k==N);UR=CqN(p(a)).*ho(q(a)).^2;up=find(UR>=0);F=UR(up);end %East-west
% if (k==-N);UR=CqmN(p(a)).*ho(q(a)).^2;up=find(UR<0);F=-UR(up);end
% if (k==1);UR=Cq1(p(a)).*ho(q(a)).^2;up=find(UR>=0);F=UR(up);end  %North-south
% if (k==-1);UR=Cqm1(p(a)).*ho(q(a)).^2;up=find(UR<0);F=-UR(up);end
% value(up)=value(up)+F*3600*24/dx./h(p(a(up)));%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;

% 
% up=[];F=0;
% if (k==N);UR=CqN(p(a)).*ho(q(a)).^2.*fTide(q(a));up=find(UR>=0);F=UR(up);end %East-west
% if (k==-N);UR=CqmN(p(a)).*ho(q(a)).^2.*fTide(q(a));up=find(UR<0);F=-UR(up);end
% if (k==1);UR=Cq1(p(a)).*ho(q(a)).^2.*fTide(q(a));up=find(UR>=0);F=UR(up);end  %North-south
% if (k==-1);UR=Cqm1(p(a)).*ho(q(a)).^2.*fTide(q(a));up=find(UR<0);F=-UR(up);end
% value(up)=value(up)+F*3600*24/dx./h(p(a(up)));%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;

% 
%BEST FEB 2025
up=[];F=0;
if (k==N);UR=CqN(p(a)).*ho(q(a)).^2;up=find(UR>=0);F=UR(up);end %East-west
if (k==-N);UR=CqmN(p(a)).*ho(q(a)).^2;up=find(UR<0);F=-UR(up);end
if (k==1);UR=Cq1(p(a)).*ho(q(a)).^2;up=find(UR>=0);F=UR(up);end  %North-south
if (k==-1);UR=Cqm1(p(a)).*ho(q(a)).^2;up=find(UR<0);F=-UR(up);end
value(up)=value(up)+F*3600*24/dx./h(p(a(up)));%.*min(fTide(p(a(up))),fTide(q(a(up))));%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;



% up=[];F=0;
% if (k==N);UR=CqN(p(a)).*ho(q(a)).^1*1;up=find(UR>=0);F=UR(up);end %East-west
% if (k==-N);UR=CqmN(p(a)).*ho(q(a)).^1*1;up=find(UR<0);F=-UR(up);end
% if (k==1);UR=Cq1(p(a)).*ho(q(a)).^1*1;up=find(UR>=0);F=UR(up);end  %North-south
% if (k==-1);UR=Cqm1(p(a)).*ho(q(a)).^1*1;up=find(UR<0);F=-UR(up);end
% value(up)=value(up)+F*3600*24/dx./h(p(a(up)));%.*min(fTide(p(a(up))),fTide(q(a(up))));%.*(fTide(p(a(up)))+fTide(q(a(up))))/2;


%value(up)=value(up)+F*3600*24/dx./h(p(a(up)))./max(0.1,fTide(p(a(up))));
%value(up)=value(up)+F*3600*24/dx./h(p(a(up))).*min(fTide(p(a(up))),fTide(q(a(up))));
end

S(p(a))=S(p(a))+value.*Fin; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value.*Fout]; %gain from the neigborh cell

end

%summary of the material that exits the cell
settling=24*3600*WS(p)./h(p).*fTide(p);%.*(h(p)>3);
settlingo=24*3600*WS(p)./h(p).*fTide(p);%.*(h(p)>3);

%sea boundary
a=find(A(p)==2);%find the co b.c.
settling(a)=0;%do not settle in the b.c. (to not change the SSC)
settlingo(a)=0;%do not settle in the b.c. (to not change the SSC)
S(p(a))=1;%to impose the b.c.

%river boundary
if computeriver==1; %if mud, handlew this as an imposed SSC
a=find(A(p)>=10 & A(p)<=19);%find the co b.c.
settling(a)=0;%do not settle in the b.c. (to not change the SSC)
S(p(a))=1;%to impose the b.c.
end

i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)+settling];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;%P=mldivide(ds2,rhs);

SSM=0*d;SSM(G>0)=full(P(G(G>0)));%rescale the matrix



%update the bed
%EmD=0*A;EmD(p)=(E(p)-SSM(p).*settling)/rbulk;
EmD=0*A;EmD(p)=(E(p)-SSM(p).*settlingo)/rbulk;














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!
p = find(A==1);%[row col]=ind2sub(size(A),p);



%Sea boundary%%%%%%%%%%%%%%%
Q=0;QR=0;
for k = [N -1 1 -N]
   
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

%only get the cell that discharge into a b.c. A==2
a=a(A(q(a))==2);

%for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;
DD=(D(p(a))+D(q(a)))/2;%.*min(fTide(p(a)),fTide(q(a)));%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
DD(isnan(DD))=0;%new Oct 2025
% if numeric==1
% DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
% elseif numeric==0
%DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
% end

%Tide
Q=Q+sum(DD.*(SSM(p(a))./h(p(a))-SSM(q(a))./h(q(a)))); %exit from that cell
%Q=Q+sum(       DD.*(  SSM(p(a))./h(p(a))./fTide(p(a))  -SSM(q(a))./h(q(a))./fTide(q(a)) )        ); %exit from that cell

% %River -FOR STANRDAD CURRENT DIRECTION
%if ebbflood==1;
if computeriver==1;
 if (k==N); QR=QR+sum( SSM(p(a))./h(p(a)).*qN(p(a)) .*(qN(p(a))>0).*(fTide(p(a))+fTide(q(a)))/2 )*sign(k);end %rigth boudnary?
 if (k==-N);QR=QR+sum( SSM(p(a))./h(p(a)).*qmN(p(a)).*(qmN(p(a))<0).*(fTide(p(a))+fTide(q(a)))/2)*sign(k);end %I think this is the left bounday
 if (k==1); QR=QR+sum( SSM(p(a))./h(p(a)).*q1(p(a)) .*(q1(p(a))>0).*(fTide(p(a))+fTide(q(a)))/2 )*sign(k);end
 if (k==-1);QR=QR+sum( SSM(p(a))./h(p(a)).*qm1(p(a)).*(qm1(p(a))<0).*(fTide(p(a))+fTide(q(a)))/2 )*sign(k);end
% if (k==N); QR=QR+sum( SSM(p(a))./h(p(a)).*qN(p(a)) .*(qN(p(a))>0))*sign(k);end %rigth boudnary?
% if (k==-N);QR=QR+sum( SSM(p(a))./h(p(a)).*qmN(p(a)).*(qmN(p(a))<0))*sign(k);end %I think this is the left bounday
% if (k==1); QR=QR+sum( SSM(p(a))./h(p(a)).*q1(p(a)) .*(q1(p(a))>0))*sign(k);end
% if (k==-1);QR=QR+sum( SSM(p(a))./h(p(a)).*qm1(p(a)).*(qm1(p(a))<0))*sign(k);end
end
%end

%River -FOR INVERTED CURRENT DIRECTION
%if ebbflood==-1;
if computeriver==1;
if (k==N); QR=QR+sum( SSM(q(a))./h(q(a)).*qmN(q(a)) .*(qmN(q(a))<0).*(fTide(p(a))+fTide(q(a)))/2)*sign(k);end
if (k==-N);QR=QR+sum( SSM(q(a))./h(q(a)).*qN(q(a)) .*(qN(q(a))>0).*(fTide(p(a))+fTide(q(a)))/2)*sign(k);end
if (k==1); QR=QR+sum( SSM(q(a))./h(q(a)).*qm1(q(a)) .*(qm1(q(a))<0).*(fTide(p(a))+fTide(q(a)))/2)*sign(k);end
if (k==-1);QR=QR+sum( SSM(q(a))./h(q(a)).*q1(q(a))  .*(q1(q(a))>0).*(fTide(p(a))+fTide(q(a)))/2)*sign(k);end
end
%end

end

QseaTide=QseaTide+dx*dt*Q/rbulk; %(note the the divided dx is for the gradient, not for the cell width!)
QseaRiver=QseaRiver+dt*3600*24*QR/rbulk;
%%%%%%%%%%%%%%%%%%%



%River mouth%%%%%%%%%%%%%%%%
for i=1:length(coRIVER)

Q=0;QR=0;
for k = [N -1 1 -N]
    
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
%exlcude the translated cell that are the b.c.
a=a(A(q(a))==10+i-1);
if (k==N | k==-N);D=Dyy;else;D=Dxx;end
DD=(D(p(a))+D(q(a)))/2;%.*min(fTide(p(a)),fTide(q(a)));%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
DD(isnan(DD))=0;%new Oct 2025
% if numeric==1%mud
% DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
% elseif numeric==0%sand
%DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
% end

%Tide
Q=Q+sum(DD.*(SSM(p(a))./h(p(a))-SSM(q(a))./h(q(a)))); %exit from that cell
%Q=Q+sum(       DD.*(  SSM(p(a))./h(p(a))./fTide(p(a))  -SSM(q(a))./h(q(a))./fTide(q(a)) )        ); %exit from that cell






%I think the way the sediment tranpostt (for River..+ Utide) is ahandles
%differnelty ad A==10 and A==2, It is not perectly symmwetric. (June 2024%comment)
% %River
%if ebbflood==1;
if computeriver==1;
if (k==N); QR=QR+sum( SSM(q(a))./h(q(a)).*qmN(q(a)) .*(qmN(q(a))>0))*sign(k);end
if (k==-N);QR=QR+sum( SSM(q(a))./h(q(a)).*qN(q(a))  .*(qN(q(a))>0))*sign(k);end
if (k==1); QR=QR+sum( SSM(q(a))./h(q(a)).*qm1(q(a)) .*(qm1(q(a))>0))*sign(k);end
if (k==-1);QR=QR+sum( SSM(q(a))./h(q(a)).*q1(q(a))  .*(q1(q(a))>0))*sign(k);end
end
%end

%if ebbflood==-1;
if computeriver==1;
if (k==N); QR=QR+sum( SSM(p(a))./h(p(a)).*qN(p(a)) .*(qN(p(a))<0))*sign(k);end
if (k==-N);QR=QR+sum( SSM(p(a))./h(p(a)).*qmN(p(a)).*(qmN(p(a))<0))*sign(k);end
if (k==1); QR=QR+sum( SSM(p(a))./h(p(a)).*q1(p(a)) .*(q1(p(a))<0))*sign(k);end
if (k==-1);QR=QR+sum( SSM(p(a))./h(p(a)).*qm1(p(a)).*(qm1(p(a))<0))*sign(k);end
end
%end






end

%pause
QmouthTide=QmouthTide-dx*dt*Q/rbulk;
QmouthRiver=QmouthRiver-dt*3600*24*QR/rbulk;
%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FLX=[QmouthRiver;QseaTide;QseaRiver;QmouthTide];

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fluxes at arbitrary cross section
if tracksedimentfluxes==1;
iii=length(fieldnames(XX));
names = fieldnames(XX);

ABASE=A;%This is the original A, used here just to remove the A==0, which are not included in the domain

for i=1:iii
eval(['XXo=XX.' names{i}  ';' ])

A=A*0;
A(XXo(:,1))=1;
A(XXo(:,2))=2;

p = find(A==1 & ABASE>0);%[row col]=ind2sub(size(A),p);


%Sea boundary%%%%%%%%%%%%%%%
Q=0;QR=0;
for k = [N -1 1 -N]
   
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

%only get the cell that discharge into a b.c. A==2
a=a(A(q(a))==2 & ABASE(q(a))>0);

%for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;
DD=(D(p(a))+D(q(a)))/2;%.*min(fTide(p(a)),fTide(q(a)));%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
% if numeric==1
% DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
% elseif numeric==0
% DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
% end

%Tide
SSC=SSM./h;%./fTide;
%SSC(SSC>10^3 & h<0.1)=NaN;
Q=Q+nansum(DD.*(SSC(p(a))-SSC(q(a)))); %exit from that cell

%River
if computeriver==1;
if (k==N); QR=QR+sum( SSM(p(a))./h(p(a)).*qN(p(a))  )*sign(k);end
if (k==-N);QR=QR+sum( SSM(p(a))./h(p(a)).*qmN(p(a)) )*sign(k);end
if (k==1); QR=QR+sum( SSM(p(a))./h(p(a)).*q1(p(a))  )*sign(k);end
if (k==-1);QR=QR+sum( SSM(p(a))./h(p(a)).*qm1(p(a)) )*sign(k);end
end

end

XXTide(i)=dx*dt*Q/rbulk; %(note the the divided dx is for the gradient, not for the cell width!)
XXRiver(i)=dt*3600*24*QR/rbulk;

end %end of loop over different cross sections XX



else%do not calcualte the fluxes at all
XXTide=NaN;
XXRiver=NaN;
end

%%%%%%%%%%%%%%%%%%%
