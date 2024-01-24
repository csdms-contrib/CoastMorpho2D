function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,B]=initializegeometry(P);

dx=20;
N=448; %x
M=157; %y

load W
W=W(1460:end-600,50:9000);
W=W(1:20:end,1:20:end)';
W=W(:,end:-1:1);


%%%%%%%%%%%cell types
A=ones(N,M);
%sea b.c.
A(end,:)=2;

%RIVER, not used in the model
rivermouthfront=[];SPCLcell=struct;SPCLcell.rivermouthfront=rivermouthfront;

%bathymetry
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
z1=-0.2-0.2-2.5+min(1112.5,+0.1+0.4-6+3.3+dx/1000*2.5*ones(N,1)*([1:M])) +min(999+8,+0.5+0.5-2.3-1+dx/1000*1.2*([1:N])'*ones(1,M));
z2=0.1+0.2+0.3+0.2-1+0.2+min(999+8,dx/1000*0.4*([1:N])'*ones(1,M));
z=min(z1,z2);  

Wg=W;
Wg(:,end)=1;%top barrier island
Wg(:,1)=1;%top barrier island
Wg(end,:)=0;%sea

zg=double(bwdist(Wg))*dx/20/3-0.5 +0.5-0.5 -0.2 -1 -0.2-0.3;%ORIGINAL
z=min(z,zg);
z=z+1.4+0.1-0.1 +0.1 +0.2+0.1+0.1+0.1 +0.1;%0.1;

%random
z=z+2*(rand(N,M)-0.5)*0.5;% *0.5   *0.2;

z(end,:)=10;

z=max(z,-1.8);%marsh cant be too high

z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%


A(W==1)=0;

msl=0;
zbedo=-z-msl;
Active=zbedo<2;%P.Trange/2;

z(A==0)=0;

B=A*0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%wave boundary conditions
%the lateral boundaries for waves
%postivie is left boundary; negative is right boudndary;
%1 is no-gradient in wave. THIS IS ALSO no-gradient in wave-indcued lateral
%sediment transport
%2 is zero wave height . Also implies no sand transport at inlet
AW=A*0;%the internal points
%right boundary
AW(:,end)=-1;
%left boundary
AW(:,1)=1;
%AW(1:50,1)=1;
%%%%%%%%%%%%%%%%%%%%%