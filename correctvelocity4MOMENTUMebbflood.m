function [Vx,Vy]=correctvelocity4MOMENTUM(N,M,periodic,h,A,MANN,dx,Ux,Uy);

MANN(isnan(MANN))=0.1;
csi=h.^(1/3)./MANN.^2/1;
D=csi.*h/dx;

%this allows to sav etime when there is a lot of puload. added may 2019
A(h<=1)=0;

p = find((A==1 | A==10));

%ADJUST FOR MOMENTUM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vx=Ux*0;Vy=Uy*0;

for k = [N -1 1 -N]  
%the translated cell
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end


if (k==1); Vx(p(a))=Vx(p(a)) -0.5*  (D(p(a))+D(q(a)))/2  /9.8.*Ux(q(a)).*Ux(p(a));end
if (k==-1);Vx(p(a))=Vx(p(a)) +0.5*  (D(p(a))+D(q(a)))/2  /9.8.*Ux(q(a)).*Ux(p(a));end
if (k==N); Vy(p(a))=Vy(p(a)) -0.5*  (D(p(a))+D(q(a)))/2  /9.8.*Uy(q(a)).*Uy(p(a));end
if (k==-N);Vy(p(a))=Vy(p(a)) +0.5*  (D(p(a))+D(q(a)))/2  /9.8.*Uy(q(a)).*Uy(p(a));end


end

Vx=sign(Vx).*sqrt(abs(Vx));
Vy=sign(Vy).*sqrt(abs(Vy));


Vx=diffuseVELOCITYxpartial(A,Vx,10,dx,h>0,0.1);
Vy=diffuseVELOCITYypartial(A,Vy,10,dx,h>0,0.1); 
