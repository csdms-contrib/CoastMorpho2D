function [z]=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,dx,dt,VEG,S,Qs,rbulk,alphaMUD);


%creepoffshore=0;

A(Active==0)=0;

%%%%%%%!!!!!!!$$$$$$
%%%%%%%!!!!!!!$$$$$$
%A(A==2)=1;%%trucco per fare creep also at the boundary


%%downslope dependnets on sediment transport
Qs=Qs/rbulk;

creep=A*0;
creep(VEG==0)=crMUD+(alphaMUD*3600*24*Qs(VEG==0));          
creep(VEG==1)=crMARSH;    

% %%%%%%%!!!!!!!$$$$$$
% %trucco per tenere il seaward basso
% creep(end-100:end,:)=creep(end-100:end,:)+100;%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=creep/(dx^2)*dt;%.*Yreduction;

%consider the pond cells as A==0
%A(S==1)=0; %DO NOT CREEP INTO PONDS!!!! NOT USED ANYMOREEEE!!!

G=0*z;
p=find(A==1);
%if creepoffshore==0;p=find(A==1);end%exclude the NOLAND CELLS(A==0)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%if creepoffshore==1;p=find(A==1 | A==2);end
%if creepoffshore==1;D(A==2)=D(A==2)*1000;end

NN=length(p);G(p)=[1:NN];rhs=z(p);[N,M]=size(G);i=[];j=[];s=[];


Spond=S;%%%%ATTENZIONE
S=0*G; %This S ia a different emanign that ponds. it is just to store values

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 
[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))==1);
%if creepoffshore==0;a=a(A(q(a))==1);end%only inclued the cells in whcih you can creep!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%if creepoffshore==1;a=a(A(q(a))==1 | A(q(a))==2);end
  
value=(D(p(a))+D(q(a)))/2; %The orginal!!!!
%value=min(D(p(a)),D(q(a))); %used ro make the wall of the marsh more vertical

value=value.*min(Yreduction(p(a)),Yreduction(q(a)));


% gradF=abs(z(p(a))-z(q(a)))/dx;
% facNL=gradF>=tan(10/180*pi);
% value=value.*facNL;


%The standrd
%DO NOT CREEP AT THE POND EDGE
value(Spond(p(a))==1 & Spond(q(a))==0) = 0;
value(Spond(p(a))==0 & Spond(q(a))==1) = 0;

% %trucco per evitare the i canali si chiudana
% %DO NOT CREEP AT THE POND EDGE
% %"it is a pond but it would be a channel, i.e., it is not a deep channel" %cit. Giulio M
% value((Spond(p(a))==1 & VEG(p(a))==1) & Spond(q(a))==0) = 0;
% value(Spond(p(a))==0 & (Spond(q(a))==1 & VEG(q(a))==1)) = 0;


S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end


% % %NOT SURE IF I NEED IT OR NOT
% % %this fix the bed elevation of the offshore condtion (A==2)
% % % if creepoffshore==1;
% % % a=find(A(p)==2);%find the co b.c.
% % % S(p(a))=0;%to impose the b.c.
% % % end

%z(A==2)
%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;z(G>0)=full(P(G(G>0)));




% %%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!
% p = find(A==2);[row col]=ind2sub(size(A),p);
% Qcreep=0;
% %Sea boundary%%%%%%%%%%%%%%%
% for k = [N -1 1 -N]
%    
% [a,q]=excludeboundarycell(k,N,M,p);
% 
% 
% %only get the cell that discharge into a b.c. A==2
% a=a(A(q(a))==1);
% %size(a)
% 
% DD=(D(p(a))+D(q(a)))/2.*min(Yreduction(p(a)),Yreduction(q(a)));
% 
% %Tide
% Qcreep=Qcreep-sum(DD.*(z(p(a))-z(q(a)))); %exit from that cell
% 
% end
% %%%%%%%%%%%%%%%%%%%
% Seacreep=Seacreep+Qcreep;

