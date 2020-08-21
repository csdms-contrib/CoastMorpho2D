function [SALT]=sedtran(d,A,SPCLcell,Dwave,DiffS,h,ho,E,ws,dx,dt,rbulk,co,Ux,Uy,FLX,fTide,Ttide,URx,URy,periodic,computeriver,computetide,kro,MUD,coMUD,Qsmouthi);


rivermouthfront=SPCLcell.rivermouthfront;


%consider the pond cells as A==0
%A(A==3)=1; %but do not update it!
%A(A==22)=2;

%A(h<=kro)=0; %eliminate the cells in which the water depth is too small
%A(rivermouthfront)=1;%the cells in front of the river mouth: let them erode if needed

p = find(A>0);%exclude the NOLAND CELLS (A==0)
G=0*d;NN=length(p);G(p)=[1:NN];
rhs=E(p)*0; %in the rhs there is already the additon of the erosion input
[N,M]=size(G);
i=[];j=[];s=[];S=0*G;

%boundary conditions imposed SSC
a=find(A==2);rhs(G(a))=co.*h(a).*fTide(a);

%iver and mud (if mud, you impose the SSC at the inlet)
if MUD==1;
a=find(A==10);rhs(G(a))=coMUD.*h(a).*fTide(a);
end

%Dturb=0.1*24*3600;  %QUESTO PER ORA COME REFERENZA, POI METTI ALTRO!!!
Dturb=5.9*0.01*h*24*3600*0;


Dxx=(Dwave*24*3600+Dturb+DiffS*Ttide/2*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, it the h in the mass balance eq.
Dyy=(Dwave*24*3600+Dturb+DiffS*Ttide/2*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*(ho>kro);%.*(hgross>0.01);
Dxx(A==10)=0;Dyy(A==10)=0; %no tidal flux at the river mouth
Dxx(rivermouthfront)=0;Dyy(rivermouthfront)=0; %no tidal flux at the river mouth
%the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]  %go over the 4 directions for the gradients

%avoid to the the cells out of the domain (risk to make it periodic...)
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

a=a(A(q(a))==1 | A(q(a))==2 | A(q(a))==10);%exlcude the translated cell that are NOLAND cells

%go over the extradiagonal direction for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;
DD=(D(p(a))+D(q(a)))/2;%.*(ho(p(a)>kro));%.*(ho(q(a))>kro); %BEST  FA UN BOTTO DI ARGINI
%DD=max(D(p(a)),D(q(a)));  %OK, NOT NECESSARY
%DD=min(D(p(a)),D(q(a))); %BAD
    
%diagonal
Fin=0*DD;Fin(A(p(a))==1)=1; % (A(q(a))==1 | A(q(a))==10)=1; to conserve the mass a the river mouth = no input
Fout=0*DD;Fout(A(q(a))==1)=1; %needed not to affect the b.c. -> Do not choose 2 and 1p
value=DD./h(p(a))./fTide(p(a));

%river flow
if computeriver==1
if (k==N);UR=URx(p(a));up=find(UR>0);F=UR(up);end
if (k==-N);UR=URx(p(a));up=find(UR<0);F=-UR(up);end
if (k==1);UR=URy(p(a));up=find(UR>0);F=UR(up);end
if (k==-1);UR=URy(p(a));up=find(UR<0);F=-UR(up);end
value(up)=value(up)+F*3600*24/dx;
end


S(p(a))=S(p(a))+value.*Fin; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value.*Fout]; %gain from the neigborh cell

end


%sea boundary
a=find(A(p)==2);%find the co b.c.
S(p(a))=1;%to impose the b.c.

%river boundary
if MUD==1; %if mud, handlew this as an imposed SSC
a=find(A(p)==10);%find the co b.c.
S(p(a))=1;%to impose the b.c.
end

i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;
%P=mldivide(ds2,rhs);

SALT=0*d;SALT(G>0)=full(P(G(G>0)));%rescale the matrix

