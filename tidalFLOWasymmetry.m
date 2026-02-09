function [U,UxNOADV,UyNOADV,U1,Um1,UN,UmN,q1,qm1,qN,qmN,Uebb,Ufld]=flowBasin(A,manning,h,ho,dx,DH,T,periodic,kro,DD1,DDN,DDUlimit,tidalnonlinearflow,Ubase,fTide,ahBNK);
%facADVh=0.5;
%factorU=factorU.^2;
%manning=manning./factorU;
%figure
%imagesc(hHR-h)
%pause

%hHR=max(hHR,0.1);
%h=max(h,0.1);

%h(h>kro)=max(h(h>kro),0.1);
%hHR(hHR>kro)=max(hHR(hHR>kro),0.1);

hHR=h;

q1=0*A;qm1=0*A;qN=0*A;qmN=0*A;

% UtideREF=0.5;
% %UtideREF=0.2;
% if tidalnonlinearflow==1
% %Uo=max(0.1,Ubase);
% 
% Uh=0.1*sqrt(9.81*h);
% %Uo=max(0.1,sqrt(Uh.*UtideREF));
% Uo=max(0.1,Uh);
% %Uo=1;
% 
% 
% %Uo=max(0.1,Uh);
% %Uo=1;
% %Uo=(Ubase+UtideREF)/2;
% %Uo=(Ubase+1)/2;
% %Uo=(Ubase+max(0.2,0.05*sqrt(9.81*h)))/2;
% else
% Uo=1;
% end
Uo=Ubase;%.*fTide;%./fTide.^0.1;
%Uo=max(Uo,0.01);

%csi=h.^(1/3)  ./(manning.^2.*Uo);%original
%Icsi=max(0.1,h).^(1/3)  ./(manning.^2.*Uo);



%hold=h;
%h=hHR;



%if boundary free slip
%h(A==0)=NaN;
%if boundary no slip
hHR(A==0)=kro;
% 
% hmin=min(hHR,[hHR(1,:); hHR(1:end-1,:)]);
% hmin=min(hmin,[hHR(2:end,:); hHR(end,:) ]);
% hmin=min(hmin,[hHR(:,1) hHR(:,1:end-1)]);
% hmin=min(hmin,[hHR(:,2:end) hHR(:,end) ]);
hmin=min(hHR,[hHR(1,:)*0+kro; hHR(1:end-1,:)]);
hmin=min(hmin,[hHR(2:end,:); hHR(end,:)*0+kro ]);
hmin=min(hmin,[hHR(:,1)*0+kro hHR(:,1:end-1)]);
hmin=min(hmin,[hHR(:,2:end) hHR(:,end)*0+kro ]);
%hCLEAN=h;

FahBNK=ahBNK;%exp(-dx*ahBNK);%0.2;  %this is large (0.8) for small dx. HOW MUCH EFFECT FROM THE SIDES. LARGE=more friction

%option1
%hF=hmin*FahBNK+h*(1-FahBNK);

% %option2
% DH=h-hmin;
 %hF=hmin+DH*(1-FahBNK);
 %hF=hmin+min(DH,max(0.1,DH*(1-FahBNK)));

% %  %option3
  hF=hmin*FahBNK+hHR*(1-FahBNK);
% hF(hF<0.2)=min(0.2,h(hF<0.2));
% hF(hF<0.1)=min(0.1,h(hF<0.1));

%figure;imagesc(h);pause
%hF=max(hF,0.1);




%Solid boundary friction
Amin=min(A,[A(1,:); A(1:end-1,:)]);
Amin=min(Amin,[A(2:end,:); A(end,:) ]);
Amin=min(Amin,[A(:,1) A(:,1:end-1)]);
Amin=min(Amin,[A(:,2:end) A(:,end) ]);
hF(Amin==0)=min(1,hF(Amin==0));



%hF=h;%for friction


%h=hold;%for continuity

Icsi=hF.*max(0.1,hF).^(1/3)  ./(manning.^2.*Uo);
%Icsi=hF.*max(1).^(1/3)  ./(manning.^2.*Uo);

%D=csi.*h.^2/(dx^2);

G=0*h;a=find(A~=0);NN=length(a);G(a)=[1:NN];
rhs=ones(NN,1).*DH(a)/(T/2*3600*24); %in m/s!!!

[N,M] = size(G);i=[];j=[];s=[];

%boundary conditions imposed water level
a=find(A==2);
i=[i;G(a)]; j=[j;G(a)]; s=[s;ones(size(a))];rhs(G(a))=0;%water level zero




S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A==1 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

% if abs(k)==1
% DDU=h.*DD1/9.81;
% elseif abs(k)==N
% DDU=h.*DDN/9.81;
% end
%DDUm=(DDU(p(a))+DDU(q(a)))/2;
%DDUm=sign(DDU(p(a))+DDU(q(a)))/2.*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);


%Icsim=(Icsi(p(a))+Icsi(q(a)))/2;%.*max(factorU(p(a)),factorU(q(a)));
Icsim=min(Icsi(p(a)),Icsi(q(a)));%
%Icsim=0.5*((Icsi(p(a))+Icsi(q(a)))/2 + min(Icsi(p(a)),Icsi(q(a))));%.*max(factorU(p(a)),factorU(q(a)));
%Icsim=max(Icsi(p(a)),Icsi(q(a)));

%hm=0.5*(min(h(p(a)),h(q(a)))  + (h(p(a))+h(q(a)))/2 );
hm=min(h(p(a)),h(q(a)));%
%hm=(h(p(a))+h(q(a)))/2;
%hm=0.5*(hm+hmin);


%val=min(0,max(-DDUlimit*Icsim,DDUm));
%DD=1./(Icsim+val)   .*hm.^2 /(dx^2);
%DD=1./Icsim   .*hm.^2 /(dx^2);%origional
DD=Icsim   .*hm /(dx^2);%new
%DD=1./(Icsim+max(-DDUlimit*Icsim,DDUm))   .*hm.^2 /(dx^2);


S(p(a))=S(p(a))+DD; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-DD]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));
P(A==2)=0;  %need when swtinching q and p





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1=0*A;Um1=0*A;UN=0*A;UmN=0*A;
p = find(A==1 | A==2 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]
if k==N; a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1; a=find(row+1<=N);end;
q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

% if abs(k)==1
% DDU=h.*DD1/9.81;
% elseif abs(k)==N
% DDU=h.*DDN/9.81;
% end
%DDUm=(DDU(p(a))+DDU(q(a)))/2;
%DDUm=sign(DDU(p(a))+DDU(q(a)))/2.*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);

%Icsim=(Icsi(p(a))+Icsi(q(a)))/2;%.*max(factorU(p(a)),factorU(q(a)));
Icsim=min(Icsi(p(a)),Icsi(q(a)));%
%Icsim=0.5*((Icsi(p(a))+Icsi(q(a)))/2 + min(Icsi(p(a)),Icsi(q(a))));
%Icsim=max(Icsi(p(a)),Icsi(q(a)));


%hm=0.5*(min(h(p(a)),h(q(a)))  + (h(p(a))+h(q(a)))/2 );
hm=min(h(p(a)),h(q(a)));%originla
%hm=(h(p(a))+h(q(a)))/2;

%hm=max(h(p(a)),h(q(a)));%

%hm=0.5*(hm+hmin);

%DD=Icsim.*hm /(dx^2) *dx  .*min(h(p(a)),h(q(a)))./max(h(p(a)),h(q(a))) ;%.*(fTide(p(a))+fTide(q(a)))/2;

%DD=Icsim.*hm /(dx^2) *dx  .*min(h(p(a)),h(q(a)))./hm ;%.*(fTide(p(a))+fTide(q(a)))/2;

%DD=Icsim /(dx^2) *dx  .*min(h(p(a)),h(q(a)));%.*(fTide(p(a))+fTide(q(a)))/2;
DD=Icsim.*hm /(dx^2) *dx  ;%.*(fTide(p(a))+fTide(q(a)))/2;


%hmean=(h(p(a))+h(q(a)))/2;
%hmean=hm;

%Icsim=max(Icsi(p(a)),Icsi(q(a)));
%Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
%DDh=Icsim /(dx^2) *dx;
hmean=(h(p(a))+h(q(a)))/2;
%hmean=hm;


if (k==1); q1(p(a))=q1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==-1); qm1(p(a))=qm1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==N); qN(p(a))=qN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;  
elseif (k==-N); qmN(p(a))=qmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
end


if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean;  
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
end


end




%Uy=(U1+Um1)/2;
%Ux=(UN+UmN)/2;





Uy=max(abs(U1),abs(Um1)).*sign(U1+Um1);%originla
Ux=max(abs(UN),abs(UmN)).*sign(UN+UmN);%originla
 
% Uy=(U1+Um1)/2;
% Ux=(UN+UmN)/2;
% Uy=(Uy + sign(Uy).*max(abs(U1),abs(Um1)))/2;
% Ux=(Ux + sign(Ux).*max(abs(UN),abs(UmN)))/2;






%Ux(A>=10 & A<=19)=(q1(A>=10 & A<=19)+qm1(A>=10 & A<=19))/2./h(A>=10 & A<=19);
%Uy(A>=10 & A<=19)=(qN(A>=10 & A<=19)+qmN(A>=10 & A<=19))/2./h(A>=10 & A<=19);
% 

%  Uy=(q1+qm1)/2; 
%  Ux=(qN+qmN)/2;
%  Uy=0.5*(sign(Uy).*max(abs(q1),abs(qm1))+Uy)./h;
%  Ux=0.5*(sign(Ux).*max(abs(qN),abs(qmN))+Ux)./h;
 


%Uy=0.5*(min(abs(U1),abs(Um1)).*sign(U1+Um1) +(U1+Um1)/2);%originla
%Ux=0.5*(min(abs(UN),abs(UmN)).*sign(UN+UmN) +(UN+UmN)/2);%originla

U=sqrt(Ux.^2+Uy.^2);
%U(A==2)=0;

%figure;imagesc(P);pause

% Amin=min(A,[A(1,:); A(1:end-1,:)]);
% Amin=min(Amin,[A(2:end,:); A(end,:) ]);
% Amin=min(Amin,[A(:,1) A(:,1:end-1)]);
% Amin=min(Amin,[A(:,2:end) A(:,end) ]);
%U(Amin==0)=0;






%655664362435367878r44trfvrrg
%655664362435367878r44trfvrrg
%655664362435367878r44trfvrrg
%655664362435367878r44trfvrrg
UxNOADV=Ux;
UyNOADV=Uy;


%UxR=[Ux(1,:); Ux(1:end-1,:)];
%UxL=[Ux(2:end,:); Ux(end,:)];
%UyU=[Uy(:,1) Uy(:,1:end-1)];
%UyD=[Uy(:,2:end) Uy(:,end)];
UxR=[Ux(:,1) Ux(:,1:end-1)];
UxL=[Ux(:,2:end) Ux(:,end)];
UyU=[Uy(1,:); Uy(1:end-1,:)];
UyD=[Uy(2:end,:); Uy(end,:)];




diry=(q1+qm1)/2;
dirx=(qN+qmN)/2;
a=find(dirx>0);
b=find(dirx<0);
%Ux(a)=UxR(a);
%Ux(b)=UxL(b);
Ux(a)=(UxR(a)+UxNOADV(a))/2;
Ux(b)=(UxL(b)+UxNOADV(b))/2;

a=find(diry>0);
b=find(diry<0);
%Uy(a)=UyU(a);
%Uy(b)=UyD(b);
Uy(a)=(UyU(a)+UyNOADV(a))/2;
Uy(b)=(UyD(b)+UyNOADV(b))/2;

Uebb=sqrt(Ux.^2+Uy.^2);
%Uebb=U;



Ux=A*0;Uy=A*0;
diry=-(q1+qm1)/2;
dirx=-(qN+qmN)/2;
a=find(dirx>0);
b=find(dirx<0);
%Ux(a)=UxR(a);
%Ux(b)=UxL(b);
Ux(a)=(UxR(a)+UxNOADV(a))/2;
Ux(b)=(UxL(b)+UxNOADV(b))/2;

a=find(diry>0);
b=find(diry<0);
%Uy(a)=UyU(a);
%Uy(b)=UyD(b);
Uy(a)=(UyU(a)+UyNOADV(a))/2;
Uy(b)=(UyD(b)+UyNOADV(b))/2;

Ufld=sqrt(Ux.^2+Uy.^2);
%Ufld=U;

%%%%%%%%%%%%%%%%%%%%%%%%%%
