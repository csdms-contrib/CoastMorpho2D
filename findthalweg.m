function thwg=findthalweg(d,dlim,Flim)
%dlim=0.1;
%Flim=50;

d(d<dlim)=NaN;
dem=100-d;[m,n]=size(d);Y=[1:m]'*ones(n,1)';X=ones(m,1)*[1:n];

DEM = GRIDobj(X,Y,dem);

FD = FLOWobj(DEM,'preprocess','fill');
A = flowacc(FD);
F=A.Z;
thwg=0*F;thwg(F>Flim)=1;

figure;imagesc(F)