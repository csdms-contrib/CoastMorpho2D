clear;close all;clc
k=0.4;
ws=0.02;
tau=1;
us=sqrt(tau/1030);
p=ws/0.4/us
%p=1/32;

zo=0.01;
D=1;
dz=0.0001;
z=[zo:dz:D];
%C=(z/zo).^-p;
C=((1-zo)./zo./((1-z)./z)).^-p;

figure;plot(C,z)

sum(C)*dz/D;

(1/zo)^-p/(-p+1)*(D^(-p+1)-zo^(-p+1))/D;
xlim([0 1])