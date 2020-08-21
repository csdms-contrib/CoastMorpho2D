function [si,zi,Ai]=getstrat2plot(zb,flyrU,flyrb,flyr,nlyr,dlyr,plyr,Y,N,M,Mi,Yb,dx,transx_y);

%to avoid overshoot in the plotting
plyr=plyr+floor(min(0,Y));
flyr(:,:,plyr+1:end)=NaN;
Y=max(Y,0);


% if transx_y==1
% Q=squeeze(isfinite(flyr(:,Mi,:)));Q=cumsum(Q,2);
% zi=[zb(:,Mi)  zb(:,Mi)+Yb(:,Mi)   zb(:,Mi)*ones(1,nlyr)+Yb(:,Mi)*ones(1,nlyr)+Q*dlyr   zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+Y(:,Mi)];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
% si=(ones(nlyr+3,1)*[0:N-1]*dx)';
% Ai=[flyrb(:,Mi)  squeeze(flyr(:,Mi,:))   flyrU(:,Mi)   flyrU(:,Mi)*0];
% 
% elseif transx_y==2;
% Q=squeeze(isfinite(flyr(Mi,:,:)));Q=cumsum(Q,2);
% zi=[0*Yb(Mi,:)'   Yb(Mi,:)'   Yb(Mi,:)'*ones(1,nlyr)+Q*dlyr   Yb(Mi,:)'+plyr(Mi,:)'*dlyr+Y(Mi,:)'];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
% si=(ones(nlyr+3,1)*[0:M-1]*dx)';
% Ai=[flyrb(Mi,:)'  squeeze(flyr(Mi,:,:))   flyrU(Mi,:)'   flyrU(Mi,:)'*0];
% 
% end

%Y(:,Mi)
%pause
%with the white on top
if transx_y==1
Q=squeeze(isfinite(flyr(:,Mi,:)));Q=cumsum(Q,2);
zi=[zb(:,Mi)  zb(:,Mi)+Yb(:,Mi)   zb(:,Mi)*ones(1,nlyr)+Yb(:,Mi)*ones(1,nlyr)+Q*dlyr   zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+Y(:,Mi) zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+max(0,Y(:,Mi))];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
si=(ones(nlyr+4,1)*[0:N-1]*dx)';
Ai=[flyrb(:,Mi)  squeeze(flyr(:,Mi,:))   flyrU(:,Mi)   flyrU(:,Mi)*0-0.01  flyrU(:,Mi)*0-0.01];

elseif transx_y==2;
Q=squeeze(isfinite(flyr(Mi,:,:)));Q=cumsum(Q,2);
% zi=[0*Yb(Mi,:)'   Yb(Mi,:)'   Yb(Mi,:)'*ones(1,nlyr)+Q*dlyr   Yb(Mi,:)'+plyr(Mi,:)'*dlyr+Y(Mi,:)'];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
% si=(ones(nlyr+3,1)*[0:M-1]*dx)';
% Ai=[flyrb(Mi,:)'  squeeze(flyr(Mi,:,:))   flyrU(Mi,:)'   flyrU(Mi,:)'*0];

zi=[zb(Mi,:)'  zb(Mi,:)'+Yb(Mi,:)'   zb(Mi,:)'*ones(1,nlyr)+Yb(Mi,:)'*ones(1,nlyr)+Q*dlyr    zb(Mi,:)'+Yb(Mi,:)'+plyr(Mi,:)'*dlyr+Y(Mi,:)'  zb(Mi,:)'+Yb(Mi,:)'+plyr(Mi,:)'*dlyr+max(0,Y(Mi,:)')];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
si=(ones(nlyr+4,1)*[0:M-1]*dx)';
Ai=[flyrb(Mi,:)'  squeeze(flyr(Mi,:,:))   flyrU(Mi,:)'   flyrU(Mi,:)'*0-0.01  flyrU(Mi,:)'*0-0.01];
end
