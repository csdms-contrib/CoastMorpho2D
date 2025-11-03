function [Umi,Tpseai,HS,F,PWTOT,QsWslope_seaTOT]=SeaWaves(A,angle,hwSea_lim,range,wind,VEG,ndir,Nhseawave,dx,z,msl,tempdeltaMSL,hpRIV,ws1,fTide,extraHseaimposed,addextrafetch,extrafetch);
rhos=2650;ss=1.5728;

z=z-tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!

[N,M]=size(A);
Lbasin=0;%1000/dx;
Fetchlim=0;%max(50,dx*2);%dx*2;%600;%dx*2*10;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PWTOT=0;
QsWslope_seaTOT=0;

for i=1:Nhseawave
    
    dHW=max(0,msl+hpRIV+range/2-z);%water depth at MHW
    
    if Nhseawave>1
        h=max(0,dHW-range*(i-1)/(Nhseawave-1));
    else
        h=max(0,dHW);
    end
    
    MASK=0*A+1;MASK(h<=hwSea_lim | A==0 | VEG==1)=0;
    %The standard way
    if addextrafetch==0
        F=calculatefetch(MASK,ndir,dx,angle);
    elseif addextrafetch==1
        F=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch,Lbasin,MASK);
    end
    Fo=F;F(Fo<=Fetchlim)=0;
    %diffuse the fetch field
    alphadiffusefetch=1;   %messo 10 for the VCR wave validation 10;%0;%%%QUESTO ERA 1 FINO AD APRILE 23 2018!!!!!
    F=diffusefetchPERPEND(MASK,F,alphadiffusefetch,dx,angle);
    F(Fo<=Fetchlim | MASK==0)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    a=find(Fo>Fetchlim & h>hwSea_lim & F>0 & MASK==1);%h>dlo & %a=find(Fo>dx*2);%h>dlo & %a=find(h>dlo);
    D=h(a);
    Ff=F(a);
    
    TP=0*h;HS=0*h;
    %windMAP=wind*(0.5+min(0.5,0.5*Ff/2000));%after 2 km it adaptes to full wind speed.
    %windMAP=wind+Ff*0;
    [Hs,Tp]=YeV(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
    %[Hs,Tp]=YeV(Ff,windMAP,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
    HS(a)=Hs;TP(a)=Tp;TP(TP==0)=1;
    
    if extraHseaimposed==1
        HS=getextraHsea(HS,h,MASK,dx,N,M);
    end
    
    HS(h<hwSea_lim)=0;
    TP(h<hwSea_lim)=2;
    
    kwave=0*h;kk=wavek(1./TP(a),h(a));kwave(a)=kk;
    kwave(kwave==0 | h<hwSea_lim)=1;
    
    Um=pi*HS./(TP.*sinh(kwave.*h));
    cg=(2*pi./kwave./TP)*0.5.*(1+2*kwave.*h./(sinh(2*kwave.*h)));
    PW=cg*1030*9.8.*HS.^2/16;
    
    Um(h<hwSea_lim)=0;
    PW(h<hwSea_lim)=0;
    
    Um(MASK==0)=0;
    PW(MASK==0)=0;
    Umi(:,:,i)=Um;
    Tpseai(:,:,i)=TP;
    
    
    PWTOT=PWTOT+PW;
    
    QsWslope_sea=WaveSedimentTransport(HS,h,kwave,rhos,N,M,TP,dx,ss,ws1,hwSea_lim,NaN);
    QsWslope_sea(HS==0)=0;
    QsWslope_seaTOT=QsWslope_seaTOT+QsWslope_sea;
    
    %figure;imagesc(HS);caxis([0 0.5]);pause
    %figure;imagesc(F);caxis([0 10000]);pause
end

PWTOT=PWTOT/Nhseawave;
QsWslope_seaTOT=QsWslope_seaTOT/Nhseawave;

