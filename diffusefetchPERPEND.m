function F=diffusefetchPERPEND(A,F,alpha,dx,angle);

%figure;imagesc(F);pause

%D=A*(alpha*3600*24)/(dx^2)*100.*(min(1,F/2000)).^2;
%D=A*(alpha*3600*24)/(dx^2).*(min(1,F/1000)).^2;
D=A*(alpha*3600*24)/(dx^2).*(min(1,F/5000)).^2;
%figure;imagesc(D);pause
%D=A*(alpha*3600*24)/(dx^2)*20.*(min(1,F/5000)).^2;
%D=A*(alpha*3600*24)/(dx^2)*25.*(min(1,F/5000)).^2;
%max(F,500)/500*10;

% D=A*(alpha*3600*24)/(dx^2).*F.^2/1000/1000;%max(F,500)/500*10;
% Dy=D*(1+999*abs(cos(angle/180*pi)));
% Dx=D*(1+999*abs(sin(angle/180*pi)));

%D=A*(alpha*3600*24)/(dx^2)*100.*F/1000;%max(F,500)/500*10;
Dy=D.*(1+9*abs(cos(angle/180*pi)));
Dx=D.*(1+9*abs(sin(angle/180*pi)));


%Dy=D.*abs(cos(angle/180*pi));
%Dx=D.*abs(sin(angle/180*pi));


G=0*F;
p=find(A==1);%exclude the NOLAND CELLS (A==0)
NN=length(p);G(p)=[1:NN];rhs=F(p);[m,n]=size(G);i=[];j=[];s=[];

S=0*G;
[row col]=ind2sub(size(A),p);
for k = [m -1 1 -m]
%avoid to the the cells out of the domain (risk to make it periodic...)
if k==m;aa=find(col+1<=n);end;if k==-m;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=m);end;

q=p+k;%the translated cell
a=aa(A(q(aa))==1 );%only inclued the cells in whcih you can creep to
% 
if abs(k)==1;value=min(Dx(p(a)),Dx(q(a)));end
if abs(k)==m;value=min(Dy(p(a)),Dy(q(a)));end
%if abs(k)==1;value=(Dx(p(a))+Dx(q(a)))/2;end
%if abs(k)==m;value=(Dy(p(a))+Dy(q(a)))/2;end

S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;F(G>0)=full(P(G(G>0)));


