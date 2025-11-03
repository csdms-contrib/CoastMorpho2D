clear;close;clc
% load PIElarge
% d=PIElarge;
% figure;imagesc(d)
% colormap('jet')
% set(gca,'YDir','normal')
% 
% load di
% d=di;
% figure;
% figure;imagesc(d)
% colormap('jet')
% set(gca,'YDir','normal')
% caxis([-1 1])
% 
% 
% [x,y] = ginput;
% X=[1:length(d(1,:))];
% Y=interp1(x,y,X,'linear','extrap');
% for i=1:length(d(1,:))
%     d(1:floor(Y(i)),i)=999;
%     %d(floor(Y(i)):end,i)=999;
% end
% 
% PIE2mBELLO=d;
% save PIE2mBELLO PIE2mBELLO

load PIE2mBELLO; d=PIE2mBELLO;
d(d==999)=NaN;
%d(299,122)=999;
%d(286,31)=NaN;
figure;imagesc(d)
colormap('jet')
set(gca,'YDir','normal')

% dx=4;%5/2;%0/1;
% %load PIElarge
% %d=PIElarge(1:2:end,1:2:end);
% %d=double(d);
% 
% load dG
% d=dG;
% 
% figure;imagesc(d);
% colormap('jet')
% set(gca,'YDir','normal')
% 
% [x,y] = ginput;
% X=[1:length(d(1,:))];
% Y=interp1(x,y,X,'linear','extrap');
% for i=1:length(d(1,:))
%     %d(1:floor(Y(i)),i)=NaN;
%     d(floor(Y(i)):end,i)=NaN;
% end
% 
% dG=d;
% save dG dG
% 
% figure;imagesc(d);
% colormap('jet')
% set(gca,'YDir','normal')
% 
% 

%     
%     
% PIEcut=d;
% save PIEcut PIEcut

% load PIEcut
% d=PIEcut;
% for i=1:length(d(1,:))
%     a=find(isfinite(d(:,i)));
%     d(1:a(1)-1,i)=999;
%     d(a(end)+1:end,i)=999;
% end
% figure;imagesc(d);
% colormap('jet')
% set(gca,'YDir','normal')
% 
% PIEcutNaN=d;
% save PIEcutNaN PIEcutNaN
