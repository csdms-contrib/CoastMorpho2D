function HS=getextraHsea(HS,h,MASK,dx,N,M);

LL=3000/dx;
HSextra=MASK*0;
%HSextra(end-100:end,:)=0.5*[0:100]'/100*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/500)*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/200)*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/500)*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/300)*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/400)*ones(1,M);
HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/400)*ones(1,M);
%HSextra(end-LL:end,:)=1*exp(-[LL:-1:0]'*dx/500)*ones(1,M);
HSextra(end,:)=0;
%HSextra(MASK==0)=0;
%figure;imagesc(HSextra');pause


%HS=HS+min(HSextra,h*0.1);
HS=HS+min(HSextra,h*0.5);