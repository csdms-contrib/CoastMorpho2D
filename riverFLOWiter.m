function [U,Ux,Uy,P,q1,qm1,qN,qmN,Utr]=flowBasin(A,manning,h,ho,dx,fTide,Q,Uo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit,ahBNK);
%facADVh=0.5;
%A(h<0.2)=0;
%h=max(0.1,h);
%Uo=A*0+1;

%h=max(h,0.1);

Q=abs(Q);
%A(A==11)=1; %the river front cells are normal cells
%to impose no cross flux in inlet

%xxx
%Uo=0.1;
%manning=0*A+n_manning;

%USE IT WHEN COMPARING TO DELFT#D HYDRO
%manning(h<0.2)=1;

A(A==22)=1;  %this behaves as normal flow %but do not update A!
%consider the pond cells as A==0
A(A==3)=1; %the isoalted pond behaves as normal cell (btu different depth...) %but do not update A!


kro=0.01;
%h(A==0)=NaN;%A cells will create no extra friction
h(A==0)=kro;%equal to kro %depth of A cells
hmin=min(h,[h(1,:); h(1:end-1,:)]);
hmin=min(hmin,[h(2:end,:); h(end,:) ]);
hmin=min(hmin,[h(:,1) h(:,1:end-1)]);
hmin=min(hmin,[h(:,2:end) h(:,end) ]);
%hCLEAN=h;
%ahBNK=0.01;%0.01;%0.2;
FahBNK=ahBNK;%exp(-dx*ahBNK);%0.2;  %this is large (0.8) for small dx. HOW MUCH EFFECT FROM THE SIDES. LARGE=more friction

%option1
%hF=hmin*FahBNK+h*(1-FahBNK);

% %option2
% DH=h-hmin;
 %hF=hmin+DH*(1-FahBNK);
 %hF=hmin+min(DH,max(0.1,DH*(1-FahBNK)));

% %  %option3
  hF=hmin*FahBNK+h*(1-FahBNK);
 %hF(hF<0.2)=min(0.2,h(hF<0.2));
% hF(hF<0.1)=min(0.1,h(hF<0.1));


% %option4
 %FahBNK=FahBNK*min(1,h/2);
 %FahBNK=exp(-dx*ahBNK./h);
 %hF=hmin.*FahBNK+h.*(1-FahBNK);


% hFLOWmin=0.1;
% hF=max(hF,hFLOWmin);



Icsi=hF.*max(0.1,hF).^(1/3)  ./(manning.^2.*Uo);
%Icsi=hF.*hF.^(1/3)  ./(manning.^2.*Uo);


G=0*A;a=find(A~=0);NN=length(a);G(a)=[1:NN];
rhs=zeros(NN,1);  
for i=1:length(Q)
    rhs(G(A==10+i-1))=Q(i)/dx;
end


[N,M] = size(G);i=[];j=[];s=[];

%boundary conditions imposed water level
a=find(A==2);
i=[i;G(a)]; j=[j;G(a)]; s=[s;ones(size(a))];rhs(G(a))=0;%water level zero

S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A==1 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells


%Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
Icsim=min(Icsi(p(a)),Icsi(q(a)));
%Icsim=max(Icsi(p(a)),Icsi(q(a)));


%hm=(h(p(a))+h(q(a)))/2;
%hm=min(h(p(a)),h(q(a)));
hm=0.5*(min(h(p(a)),h(q(a)))  + (h(p(a))+h(q(a)))/2 );

DD=Icsim.*hm /(dx^2);


if imposenocrossflux==1;
for iii=1:length(Q)  
    if abs(directionQ(iii))==1 & abs(k)==N;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    elseif abs(directionQ(iii))==2 & abs(k)==1;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    end
end
end


S(p(a))=S(p(a))+DD; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-DD]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];

ds2 = sparse(i,j,s);%solve the matrix inversion

p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));
P(A==2 | A==21)=0;  %need when swtinching q and p











%calculate flow
U1=0*A;Um1=0*A;UN=0*A;UmN=0*A;
q1=0*A;qm1=0*A;qN=0*A;qmN=0*A;
Ux=0*A;Uy=0*A;

p = find(A==1 | A==2 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]
if k==N; a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1; a=find(row+1<=N);end;
q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells


%Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
Icsim=min(Icsi(p(a)),Icsi(q(a)));
%Icsim=max(Icsi(p(a)),Icsi(q(a)));

%hm=(h(p(a))+h(q(a)))/2;
%hm=min(h(p(a)),h(q(a)));
hm=0.5*(min(h(p(a)),h(q(a)))  + (h(p(a))+h(q(a)))/2 );

DD=Icsim.*hm /(dx^2) *dx;


hmean=(h(p(a))+h(q(a)))/2;
%DDh=Icsim /(dx^2) *dx;

if imposenocrossflux==1;
for iii=1:length(Q)  
    if abs(directionQ(iii))==1 & abs(k)==N;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    elseif abs(directionQ(iii))==2 & abs(k)==1;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    end
end
end

if (k==1); q1(p(a))=q1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==-1); qm1(p(a))=qm1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==N); qN(p(a))=qN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;  
elseif (k==-N); qmN(p(a))=qmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
end

if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean;  
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hmean; 
end


end







for i=1:length(Q)  
    
if directionQ(i)==1 %up
qm1(A==10+i-1)=q1(A==10+i-1);

elseif directionQ(i)==2 %left
qmN(A==10+i-1)=qN(A==10+i-1);

elseif directionQ(i)==-2 %right
qN(A==10+i-1)=qmN(A==10+i-1);
end

end





%%%ONLY THIS WORKS FOR THE RIVER!!!!
%%option 1
 %Ux=(q1+qm1)/2./h; 
 %Uy=(qN+qmN)/2./h;
 %Ux=min(q1,qm1)./h; 
 %Uy=min(qN,qmN)/2./h;
 %Ux=max(abs(q1),abs(qm1))./h;
 %Uy=max(abs(qN),abs(qmN))./h;
 
%option 2
%Ux=(U1+Um1)/2;
%Uy=(UN+UmN)/2;



%option BBB-more incidpesd braided river?, more bank erosion?
Ux=max(abs(U1),abs(Um1));
Uy=max(abs(UN),abs(UmN));


%option  clasic smmoth
% Ux=(U1+Um1)/2;
% Uy=(UN+UmN)/2;
% Ux=(Ux + sign(Ux).*max(abs(U1),abs(Um1)))/2;
% Uy=(Uy + sign(Uy).*max(abs(UN),abs(UmN)))/2;

Ux(A>=10 & A<=19)=(q1(A>=10 & A<=19)+qm1(A>=10 & A<=19))/2./h(A>=10 & A<=19);
Uy(A>=10 & A<=19)=(qN(A>=10 & A<=19)+qmN(A>=10 & A<=19))/2./h(A>=10 & A<=19);
%  
% %UG=sqrt(Ux.^2+Uy.^2);
%  
%  


%Ux=min(abs(q1),abs(qm1))./h; %THIS IS TERRIBLE
%Uy=min(abs(qN),abs(qmN))./h; %THIS IS TERRIBLE
%Ux=max(abs(q1),abs(qm1))./h; THIS IS TERRIBLE
%Uy=max(abs(qN),abs(qmN))./h;THIS IS TERRIBLE

% 
% 
% %option 2
 %Ux=(U1+Um1)/2;
 %Uy=(UN+UmN)/2;
  %Ux=max(abs(U1),abs(Um1));
  %Uy=max(abs(UN),abs(UmN));
% 
%  
% % UxT=(U1+Um1)/2;
% % UyT=(UN+UmN)/2;
% % Ux(UxT>=0)=U1(UxT>=0);
% % Ux(UxT<0)=Um1(UxT<0);
% % Uy(UyT>=0)=UN(UyT>=0);
% % Uy(UyT<0)=UmN(UyT<0);
%  
% Ux(A>=10 & A<=19)=(q1(A>=10 & A<=19)+qm1(A>=10 & A<=19))/2./h(A>=10 & A<=19);
% Uy(A>=10 & A<=19)=(qN(A>=10 & A<=19)+qmN(A>=10 & A<=19))/2./h(A>=10 & A<=19);

% % %option 3
% % Ux=(U1+Um1)/2;
% % Uy=(UN+UmN)/2;
% % Ux(A>=10 & A<=19)=(q1(A>=10 & A<=19)+qm1(A>=10 & A<=19))/2./h(A>=10 & A<=19);
% % Uy(A>=10 & A<=19)=(qN(A>=10 & A<=19)+qmN(A>=10 & A<=19))/2./h(A>=10 & A<=19);
% % Ux=(Ux+(q1+qm1)/2./h)/2; 
% % Uy=(Uy+(qN+qmN)/2./h)/2;
% 



%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%for transport
Utr=sqrt(Ux.^2+Uy.^2);




%for fricton
 %Ux=max(abs(q1),abs(qm1))./h;
 %Uy=max(abs(qN),abs(qmN))./h;
 Ux=(q1+qm1)/2;%max(h,0.1);
 Uy=(qN+qmN)/2;%max(h,0.1);
U=sqrt(Ux.^2+Uy.^2)./h;%max(h,0.1);



% %VERY GOOD 11 12 2024
% % %Ux=(q1+qm1)/2./h; 
% % %Uy=(qN+qmN)/2./h;
%  Ux=(U1+Um1)/2;
%  Uy=(UN+UmN)/2;
% % %Ux=min(abs(q1),abs(qm1))./h;
% % %Uy=min(abs(qN),abs(qmN))./h;
%  Utr=sqrt(Ux.^2+Uy.^2);
% % 
%  Utr(A>2)=U(A>2);
%  
 
 
% %TEST 11 13 2024!!!
% UxT=(U1+Um1)/2;
% UyT=(UN+UmN)/2;
% Ux=0*A;Uy=0*A;
% Ux(UxT>=0)=Um1(UxT>=0);
% Ux(UxT<0)=U1(UxT<0);
% Uy(UyT>=0)=UmN(UyT>=0);
% Uy(UyT<0)=UN(UyT<0);
% Utr=sqrt(Ux.^2+Uy.^2);
% 
% Utr(A>2)=U(A>2);
%  
 
 
 
 
 
 
 
 
 
 


% Uo=min(U,FcrUR*sqrt(9.81*h));
% 
% Ux=Uo.*Ux./U;Ux(U==0)=0;
% Uy=Uo.*Uy./U;Uy(U==0)=0;
% U=sqrt(Ux.^2+Uy.^2);






















% 
% %this will impose more uniform discharge over the oulet
% for i=1:length(Q)
%     
%     
%     if directionQ(i)==1 %up
%     Uy(A==10+i-1)=Q(i)./h(A==10+i-1); Ux(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     qN(A==10+i-1)=0;qmN(A==10+i-1)=0;%mouth cell
%     %q1(A==10+i-1)=Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qm1(q)=Q(i);end %cell in front of the mouth
% %     
%     elseif directionQ(i)==2 %left
%     Ux(A==10+i-1)=Q(i)./h(A==10+i-1); Uy(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     q1(A==10+i-1)=0;qm1(A==10+i-1)=0;%mouth cell
%     %qN(A==10+i-1)=Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qmN(q)=Q(i);end %cell in front of the mouth
%     
%     elseif directionQ(i)==-2 %right
%     Ux(A==10+i-1)=Q(i)./h(A==10+i-1); Uy(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     q1(A==10+i-1)=0;qm1(A==10+i-1)=0;%mouth cell
%     %qmN(A==10+i-1)=-Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qN(q)=-Q(i);end %cell in front of the mouth
%    
%     end
% end

% figure
% subplot(2,1,1);
% imagesc(q1);caxis([0 82]);colormap('jet')
% subplot(2,1,2);
% imagesc(qm1);caxis([0 82]);colormap('jet')
% pause




%q=U.*h;












% q1=A*0+0.1;
% qm1=A*0+0.1;
% U=A*0+0.1;

%figure;imagesc(qm1);pause
% 
% figure;imagesc(P);
% 
% Dx=2*P-[P(1,:); P(1:end-1,:)]-[P(2:end,:); P(end,:)];
% Dy=2*P-[P(:,1) P(:,1:end-1)]-[P(:,2:end) P(:,end)];
% 
% 
% figure;
% subplot(3,1,1);imagesc(Dx);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% subplot(3,1,2);imagesc(Dy);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% subplot(3,1,3);imagesc(Dx+Dy);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% 
% 
% 




% 
%%%check for divergence free

% Dx=[U1(1,:); U1(1:end-1,:)]-[Um1(2:end,:); Um1(end,:)];
% Dy=[UN(:,1) UN(:,1:end-1)]-[UmN(:,2:end) UmN(:,end)];
% 
% figure;
% subplot(3,1,1);imagesc(Dx);caxis([-1 1]/5);colormap('jet')%;xlim([0 10]);ylim([80 120])
% subplot(3,1,2);imagesc(-Dy);caxis([-1 1]/5);colormap('jet')
% subplot(3,1,3);imagesc(Dx+Dy);caxis([-1 1]/5);colormap('jet')
% pause







% % 
% % % P1=[P(:,2:end) P(:,end)];P2=[P(:,1) P(:,1:end-1) ];P3=[P(2:end,:); P(end,:)];P4=[P(1,:); P(1:end-1,:)];
% % % D1=[D(:,2:end) D(:,end)];D2=[D(:,1) D(:,1:end-1) ];D3=[D(2:end,:); D(end,:)];D4=[D(1,:); D(1:end-1,:)];
% % % A1=[A(:,2:end) A(:,end)];A2=[A(:,1) A(:,1:end-1) ];A3=[A(2:end,:); A(end,:)];A4=[A(1,:); A(1:end-1,:)];
% % % 
% % % pp1=P(A1==0);pp2=P(A2==0);pp3=P(A3==0);pp4=P(A4==0);
% % % P1(A1==0)=pp1;P2(A2==0)=pp2;P3(A3==0)=pp3;P4(A4==0)=pp4;
% % % 
% % % pp1=D(A1==0);pp2=D(A2==0);pp3=D(A3==0);pp4=D(A4==0);
% % % D1(A1==0)=pp1;D2(A2==0)=pp2;D3(A3==0)=pp3;D4(A4==0)=pp4;
% % % 
% % % DD1=(D+D1)/2;DD2=(D+D2)/2;DD3=(D+D3)/2;DD4=(D+D4)/2;
% % % 
% % %   Ux=0.5*((P-P1).*DD1 + (P2-P).*DD2)./h*dx;
% % %   Uy=0.5*((P-P3).*DD3 + (P4-P).*DD4)./h*dx;
% % %   Ux(A==2 | A==10)=2*Ux(A==2 | A==10);
% % %   Uy(A==2 | A==10)=2*Uy(A==2 | A==10);
% % % %G=cat(3,(P-P1).*DD1,(P2-P).*DD2);Ux=nanmean(G,3);
% % % %G=cat(3,(P-P3).*DD3,(P4-P).*DD4);Uy=nanmean(G,3);
% % 


