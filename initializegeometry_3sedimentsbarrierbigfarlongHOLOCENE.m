function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry(P);

%geometric parameters
dx=200/1;
%N=140*1+200+300; %x
N=140*1+120; %x
N=N*1.4+8;
M=200;%74;%100*2;%70;
%150+50; %y

%%%%%%%%%%%cell types
A=ones(N,M);
%sea b.c.
A(end,:)=2;
%river b.c.
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


%slope=0.3/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-5-10-25+30;%sloping   0.5

%THAD GOOD ONE
%slope=0.7/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-5-10-25-50;%sloping   0.5
slope=0.7/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-5-10-25-12;%sloping   0.5
%slope=0.8/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-5-10-25-20-20-7-13;%sloping   0.5


%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping
%zin=5;z=zin*ones(N,M);%flat

slopebreak=70;
% z=[1:N];
% slope=1/1000;z(1:N-slopebreak)=[1:N-slopebreak]*slope*dx;
%slope=1/200;
%z(N-slopebreak+1:N,:)=50;
%z(N-slopebreak:N,:)=(z(N-slopebreak,1)+[1:slopebreak+1])'*ones(1,M)*slope*dx;
% z=z'*ones(1,M);%sloping

%barriers
pos=80*dx/dx;bwidth=4000/dx;  %6000
%z(end-pos:end-pos+bwidth,1:M/2-2)=-2.1;
%z(end-pos:end-pos+bwidth,M/2+2:end)=-2.1;
%z(end-pos:end-pos+bwidth,1:M/2-1)=-5.1;
%z(end-pos:end-pos+bwidth,M/2+1:end)=-5.1;


%THIS IS THE GOOD ONE
%z(end-pos-bwidth:end-pos,:)=1;%3;
%z(end-pos+1:end,:)=z(end-pos+1:end,:)+0;

%z(end-20:end,M/2-2:M/2+2)=3;

%random
randz=2*(rand(N,M)-0.5)*0.5;
z=z+randz;




%river
%mouthW=3;
%z(A==10)=P.hmouth;
%z(1:160,M/2-mouthW:M/2+mouthW)=P.hmouth;

%sea
%z(end,:)=50;%%%50;  cinquanta original
z(end,:)=50;%%%50;  cinquanta original
% z(end-1,:)=30;%%%50;  cinquanta original
%  z(end-2,:)=20;%%%50;  cinquanta original
%  z(end-3,:)=10;%%%50;  cinquanta original
%  z(end-4,:)=5;%%%50;  cinquanta original
% z(end-5,:)=50;%%%50;  cinquanta original
% z(end-6,:)=40;%%%50;  cinquanta original
% z(end-7,:)=30;%%%50;  cinquanta original
% z(end-8,:)=20;%%%50;  cinquanta original
% z(end-9,:)=10;%%%50;  cinquanta original
% %z(end-10,:)=0;%%%50;  cinquanta original
% %z(end-11,:)=-1;%%%50;  cinquanta original
z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%



%%%%MITSUE STEFANNI
% load zbedo;
% z=zbedo(1:400,:);
% z(end,:)=500;
% [N,M]=size(z);
% A=z*0+1;
% A(end,:)=2;
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;



%%TO CUT THE LANDWARD STUFF
%A=A(500:end,:);z=z(500:end,:);
randz=randz(200:end,:);
[N,M]=size(z);
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
%%%%%%%%%%%%%%%%%%%



% [N,M]=size(z);
% dx=100;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% A=interp2([0:N-1],[0:M-1]',A',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


% 
%River
% mouthW=1; %2
% A(1:10,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:10,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% SPC=A*0;SPC(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(SPC==1);clear SPC;
% SPCLcell=struct;
% SPCLcell.rivermouthfront=rivermouthfront;
% %river
% z(A==10)=P.hmouth;
% z(1:152,M/2-mouthW:M/2+mouthW)=P.hmouth;
% 
% 


% load zbedo;
% z=zbedo;
% [N,M]=size(z);
% dx=100;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
% %%%%%%%%%%%cell types
% A=ones(N,M);
% %sea b.c.
% A(end,:)=2;
% %river b.c.
% rivermouthfront=[];

%questo Tp10 Luglio 5 dopo barca
%load zR20lowres;z=-zR20lowres;
%load zR20wholet550;z=-zR20wholet550;
%load zR20wholet4000;z=-zR20wholet4000;



% %latest
% load zR20edget1250;z=-zR20edget1250;
% %add land
% G=-3/2-0.0005*dx*[50:-1:0]'*ones(1,M)+2*(rand(50+1,M)-0.5)*0.5-6-1.258-0.1;
% z=[G; z];
% [N,M]=size(z);
% A=z*0+1;
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;
% 




% %BELLO
% load zR20edget1250_R13650;z=-zR20edget1250_R13650;
% [N,M]=size(z);
% A=z*0+1;
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;




%load zR20lowrest600;z=-zR20lowrest600;
%load zR20lowresQw03t550;z=-zR20lowresQw03t550;



% %cut
% A=A(1:end-300,:);
% z=z(1:end-300,:);
% [N,M]=size(z);
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;

% %cut
% A=A(150:end,:);
% z=z(150:end,:);
% [N,M]=size(z);
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;


msl=0;
zbedo=-z-msl;
Active=zbedo<P.Trange/2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eend=N;%-10;%130;
% flyr1(1:eend,:,2:2:P.levo)=0.55;flyr2(1:eend,:,2:2:P.levo)=0.45;
% flyr1(1:eend,:,1:2:P.levo)=0.45;flyr2(1:eend,:,1:2:P.levo)=0.55;
% flyrb1(1:eend,:)=0.55;flyrb2(1:eend,:)=0.45;
% Y2(1:eend,:)=Y1(1:eend,:)/2;
% Y1(1:eend,:)=Y2(1:eend,:);
% 
% f1=Y1./(Y1+Y2);
% f2=Y2./(Y1+Y2);
% Y1=Y1+randz.*f1;
% Y2=Y2+randz.*f2;
% zb=zb-randz;

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

msl=0;