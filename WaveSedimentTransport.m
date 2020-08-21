function [QsWslope,QsWon]=wavetransport(Hs,hwave,kwave,rhos,N,M,Tp,dx,ss,ws,hwavelim,fTide)

%Hs=Hs*sqrt(2);


g=9.81;
K=rhos*3600*24*(16*0.01*0.01/(15*pi*ss*g)); 
%the rhos*3600*24 is needed to give kg/m/day instead of m3/m/s

QsWslope=zeros(N,M);QsWon=zeros(N,M);
a=find(hwave>hwavelim & kwave>0 & Tp>0);

Wtrans=K/ws*( pi*Hs(a)/sqrt(2)./(Tp(a).*sinh(kwave(a).*hwave(a))) ).^5;%Wtrans=K/ws*Uwave(a).^5 ;

%Wtrans=Wtrans.*fTide(a);


%Downslope transport
QsWslope(a)=Wtrans/ws;

%
redSTRtubr=1;%reduction in streaming velocity in turbulent boundary layer
%Along wave trasnport
QsWon(a)=Wtrans.*3.*( redSTRtubr*5*Tp(a)./(4*(2*pi./kwave(a)))  +3*Tp(a)./(4*(2*pi./kwave(a)).*(sinh(kwave(a).*hwave(a))).^2) );



















% figure;imagesc(Tp)
% pause

% %shift the inputs: d and kwave
% angledir=sign(angleswell);
% tgangle=tan(abs(angleswell)/180*pi);
% for i=1:N
%     s=floor(tgangle*i);
%     QsWon(i,:)=circshift(QsWon(i,:),[0 s*angledir]);
%     QsWslope(i,:)=circshift(QsWslope(i,:),[0 s*angledir]);
% end
% 
% %QsWon=[QsWon(:,2:end)  QsWon(:,end)];%leftorright
% %QsWon=[QsWon(:,1)  QsWon(:,1:end-1)];leftorright
% %QsWon=(QsWon+[QsWon(1,:);  QsWon(1:end-1,:)])/2;%makes it isntabilty
% QsWon=(QsWon+[QsWon(1,:);  QsWon(1:end-1,:)])/2;%makes it isntabilty
% %QsWon=([QsWon(1,:);  QsWon(1:end-1,:)]+[QsWon(2:end,:);  QsWon(end,:)])/2;%makes it flat
% 
% QsWslope=([QsWslope(1,:);  QsWslope(1:end-1,:)]+[QsWslope(2:end,:);  QsWslope(end,:)])/2;%makes it flat
% %QsWslope=([QsWslope(1,:);  QsWslope(1:end-1,:)]);%makes it isntabilty
% 
% for i=1:N
%     s=floor(tgangle*i);
%     QsWon(i,:)=circshift(QsWon(i,:),[0 -s*angledir]);
%     QsWslope(i,:)=circshift(QsWslope(i,:),[0 -s*angledir]);
% end