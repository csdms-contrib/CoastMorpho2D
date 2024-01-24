function [U,Ux,Uy,U1,Um1,UN,UmN]=flowBasin(A,manning,h,ho,dx,DH,T,periodic,kro,DD1,DDN,DDUlimit,tidalnonlinearflow,Ubase);
facADVh=0.5;

UtideREF=1;
if tidalnonlinearflow==1
Uo=(Ubase+UtideREF)/2;
%Uo=(Ubase+1)/2;
%Uo=(Ubase+max(0.2,0.05*sqrt(9.81*h)))/2;
else
Uo=1;
end

csi=h.^(1/3)  ./(manning.^2.*Uo);
Icsi=1./csi;

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

if abs(k)==1
DDU=h.*DD1/9.81;
elseif abs(k)==N
DDU=h.*DDN/9.81;
end
%DDUm=(DDU(p(a))+DDU(q(a)))/2;
DDUm=sign(DDU(p(a))+DDU(q(a)))/2.*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);
Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
hm=(h(p(a))+h(q(a)))/2;

val=min(0,max(-DDUlimit*Icsim,DDUm));
%DD=1./(Icsim+val)   .*hm.^2 /(dx^2);
DD=1./Icsim   .*hm.^2 /(dx^2);
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

if abs(k)==1
DDU=h.*DD1/9.81;
elseif abs(k)==N
DDU=h.*DDN/9.81;
end
%DDUm=(DDU(p(a))+DDU(q(a)))/2;
DDUm=sign(DDU(p(a))+DDU(q(a)))/2.*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);
Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
hm=min(h(p(a)),h(q(a)));


val=min(0,max(-DDUlimit*Icsim,DDUm));
%DD=1./(Icsim+val)   .*hm /(dx^2) *dx;
DD=1./Icsim   .*hm /(dx^2) *dx;


if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;  
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
end

end

Uy=max(abs(U1),abs(Um1)).*sign(U1+Um1);
Ux=max(abs(UN),abs(UmN)).*sign(UN+UmN);
U=sqrt(Ux.^2+Uy.^2);




