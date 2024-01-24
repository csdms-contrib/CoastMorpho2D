z(z>0.5)=1;
z(z<=0.5)=0;
%EXTRAEROSIONINCURVATURE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z1=A*0;Zm1=A*0;ZN=A*0;ZmN=A*0;

p = find(A==1 | A==10);[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%the translated cells
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

dz=z(p(a))-z(q(a));

if k==1;Z1(p(a))=dz;
elseif k==-1;Zm1(p(a))=dz;
elseif k==N;ZN(p(a))=dz;
elseif k==-N;ZmN(p(a))=dz;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(2,2,1);imagesc(Z1)
% subplot(2,2,2);imagesc(Zm1)
% subplot(2,2,3);imagesc(ZN)
% subplot(2,2,4);imagesc(ZmN)
% pause

zdc=0.3;
a=find((Z1>zdc & (ZN>zdc | ZmN>zdc))  |  (Zm1>zdc & (ZN>zdc | ZmN>zdc)) );
Ecurv=A*0;
Ecurv(a)=1;
figure;
subplot(1,2,1);imagesc(Ecurv);
subplot(1,2,2);imagesc(z);
pause