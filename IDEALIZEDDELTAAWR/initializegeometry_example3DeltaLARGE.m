function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,B,directionQ]=initializegeometry(P);

load initialbedlevel

% %geometric parameters
dx=25;
N=227;%;%110*1; %x
M=300;%200*1; %y

% 
%z=0.000375*dx*[0:N-1]'*ones(1,M)+2.5-0.5;
z=-data.Val';z=z(:,2:end-1);


[N,M]=size(z);
A=ones(N,M);

A(end,:)=2;

rivermouthfront=[];
% mouthW=1; %2
% A(1:2,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:2,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% S=A*0;S(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(S==1);clear S;

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;

%rn=2*(rand(N,M)-0.5)*0.1;
%z=z+rn;

z(A==0)=NaN;
z(A==2)=100;

%%%%%%%%%%%%%%%%%%

% 
% %%nearshoremask
% dM=120*250/dx;
% for i=M-dM:M
%     A(end-min(15000/dx,floor(0.75*(dM-(M-i)))):end,i)=2;
% end
% 
% %%New Orelans
% dM=170*250/dx;
% for i=M-dM:M
%     A(1:min(14000/dx,floor(1*(dM-(M-i)))),i)=0;
% end
% 
% %%Deltafarms
% dN=53000/dx;
% for i=32500/dx:dN
%     A(i,1:min(104000/dx,floor(0.5*(dN-(N-i)))))=0;
% end
% 
% %%Marina Near mid diversion
% for i=13500/dx:24000/dx
%     A(i,M-25500/dx:M)=0;
% end
% dN=33000/dx;
% for i=24000/dx:dN+1000/dx
%     A(i,end+12000/dx+(floor(2*(dN-(N-i)))):end)=0;
% end



%A(:,1)=0;
%A(:,end)=0;
A(1,:)=0;
% 
% %Left side
% A(M/2-1:M/2+1,1)=10;
% z(A==10)=P.hmouth(1);
% %z(1:5,N/2-1:N/2+1)=P.hmouth(1);

%north side
A(1,M/2-4:M/2+4+1)=10;

A(1:21,1:-1+M/2-4)=0;
A(1:21,M/2+4+1+1:end)=0;
%A(1:5,1:M/2-1-1)=0;
%A(1:5,M/2+1+1:end)=0;
z(A==10)=P.hmouth(1);
%z(2:10+20,M/2-4+3:M/2+4+1-3)=P.hmouth(1)+3;
%%z(2:10+20,M/2-4+4:M/2+4+1-4)=P.hmouth(1)+6;

% %Right side
% A(M/2-1:M/2+1,end)=11;
% z(A==11)=P.hmouth(2);
% %z(1:5,N/2-1:N/2+1)=P.hmouth(1);


%lateral open boundary
%A(22:end,1)=2;
%A(22:end,end)=2;


%figure;imagesc(A);pause
%z=z*0+P.hmouth(1);

msl=0;
zbedo=-z-msl;
Active=zbedo<999;


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
