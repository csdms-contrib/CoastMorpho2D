function Qs=QsTOT(Ueff,ss,g,d50,rhos,h,Ds)

Me=Ueff/sqrt(ss*g.*d50);
Qsb=0.015*rhos.*h.*(d50./h).^(1.2).*Me.^1.5; %kg/s/m
Qss=0.008*rhos*d50.*Me.^2.4.*Ds.^-0.6; %kg/s/m;


%if using the thrhold for suspension
%Qss=0.03*rhos*d50.*Me.^2.*Ds.^-0.6; %kg/s/m;

Qs=(Qsb+Qss);
