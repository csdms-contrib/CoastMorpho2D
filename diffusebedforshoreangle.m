function z=bedcreepponds(z,A,MASK,cr,dx,periodic);

creep=MASK*cr;

D=creep*3600*24/(dx^2);%.*Yreduction;


%consider the pond cells as A==0
%A(S==1)=0; %DO NOT CREEP INTO PONDS!!!!

G=0*z;
p=find(A==1);%exclude the NOLAND CELLS (A==0)
NN=length(p);G(p)=[1:NN];rhs=z(p);[N,M]=size(G);i=[];j=[];s=[];


S=0*G; %This S ia a different emanign that ponds. it is just to store values
[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 
%avoid to the the cells out of the domain (risk to make it periodic...)
%if k==N;aa=find(col+1<=M);end;if k==-N;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=N);end;
%q=p+k;%the translated cell
%a=aa(A(q(aa))==1);%only inclued the cells in whcih you can creep to

    %boundary cells
    if periodic==0
    [a,q]=excludeboundarycell(k,N,M,p);%NOGRADIENT
    elseif periodic==1;
    [a,q]=periodicY(k,N,M,p); %for the long-shore
    end
    a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells | A(q(a))==2

   


value=(D(p(a))+D(q(a)))/2;%.*facNL;



S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;z(G>0)=full(P(G(G>0)));


