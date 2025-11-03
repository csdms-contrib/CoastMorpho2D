function [bio]=plantbiomass_lateralexpansion(z,S,A,N,M,dx,dt,plantexpansionrate,bio,Active,Zlev,dBseed);

p=find(bio>0);%find the existing ponds

[row col]=ind2sub([N M],p);

for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;
if k==-N;a=find(col-1>0);end;
if k==-1;a=find(row-1>0);end;
if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
%er=a(A(q(a))==1 & S(q(a))==0 ); %find the neightr cell that is a non-pond cell,  %orignal AWR
%er=a(A(q(a))==1 & bio(q(a))==0  & Active(q(a))==1); %find the neightr cell that is a non-pond cell, and is active
er=a(A(q(a))==1 & bio(q(a))==0  & Active(q(a))==1  & Zlev(q(a))>dBseed); %Modified Spet 4 2025



%rng('shuffle');
r=rand(length(er),1);% you only need rng at the beginnign of the loop
aa=find(r<plantexpansionrate*dt/365/dx);

%S(q(er(aa)))=1;


%you go over all pond cells (p) and you allow them to expand along all (up to
%four) edges with a marsh cell (q). So the cell that eventally erodes is q (not p)
for i=1:length(aa)

    bio(q(er(aa(i))))=0.1;%also update the pond ID

end
end