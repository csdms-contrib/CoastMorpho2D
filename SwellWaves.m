function [H,Uwave]=wave2Dsparse(A,AW,Ho,N,M,d,T,kwave,dx,periodic,angle,gridDIR);

Holateral=Ho;
%OCIO YOU IMPOSED PERIODIC==1 in chsoign the q+k later on (foe Delaware
%trick)

alpha=0;%10;%2.5; %adimensional
stretchdx=1/cos(angle/180*pi);%when you stretch the coordinates, you also stretch the dx!!!

d(d<0 | A==0)=0;
kwave(d<=0.1)=1;

%Flip the domain upside down depending on the orientation of the grid   
if gridDIR==-1; 
    A=A(end:-1:1,:);
    AW=AW(end:-1:1,:);
    d=d(end:-1:1,:);
    kwave=kwave(end:-1:1,:);
end
    
%Flip the angle direction depending on the orientation of the grid  
angle0=angle;
angle=angle*gridDIR;


%Shear the grid inputs: d and kwave
angledir=-sign(angle)*gridDIR;
tgangle=tan(abs(angle)/180*pi);
for i=1:N
    s=floor(tgangle*i);
    d(i,:)=circshift(d(i,:),[0 s*angledir]);
    kwave(i,:)=circshift(kwave(i,:),[0 s*angledir]);
    AW(i,:)=circshift(AW(i,:),[0 s*angledir]);
end

%Repad the lateral boundariesif periodic==0
%DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
if periodic==0  
    for i=1:N-1    
    
    if angle0>0 %wave right to left. Use the right boundary
        a1=find(AW(i+1,:)==-1 | AW(i+1,:)==-2);
        a2=find(AW(i,:)==-1 | AW(i,:)==-2);
        if length(a2)>0  %this line introduced to put AW not along all side boundary, only a piece and allow some AW=0 along sides
        if a2>=a1;AW(i,a1:a2)=AW(i,a2(1));else;AW(i,1:a2)=AW(i,a2(1));AW(i,a1:M)=AW(i,a2(1));end
        end
        
    elseif angle0<0 %wave left to right. Use the left boundary
        a2=find(AW(i+1,:)==1 | AW(i+1,:)==2);
        a1=find(AW(i,:)==1 | AW(i,:)==2);
        if length(a1)>0
        if a2>=a1;AW(i,a1:a2)=AW(i,a1(1));else;AW(i,1:a2)=AW(i,a1(1));AW(i,a1:M)=AW(i,a1(1));end
        end
    end
    
    end
end



H=zeros(N,M);
H(end,:)=Ho;

sigma=2*pi/T;
c=(2*pi/T)./kwave;
cg=0.5*sigma./kwave.*(1+2*kwave.*d./sinh(2*kwave.*d));cg(d<=0.1)=1;
    
%bed friction
bedfr=1*0.038*(sigma./(9.81*sinh(kwave.*d))).^2;bedfr(d<=0.1)=1;
%bedfr=0.067*(sigma./(9.81*sinh(kwave.*d))).^2;bedfr(d<=0.1)=1;

%wave breaking crietion
%Hcr=0.14*2*pi./kwave.*tanh(kwave.*d);Hcr(d<=0.1)=0;
Hcr=0.73*d;Hcr(d<=0.1)=0.01;


%The initial E value at the sea boundary
E=H(end,:).^2/16;


    
    
p = [1:M]';%exclude the NOLAND CELLS (A==0)

%%%%%%%%%%%%%%%%%%INITIALIZE THE INDEXES AND VECTOES
ii=[];jj=[];
for k=[-1 1];
    
%if periodic==0;     q=p+k;if k==1;a=[1:M-1];elseif k==-1;a=[2:M];end 
%elseif periodic==1; q=1+mod(p+k-1,M);a=[1:M];
%end
q=1+mod(p+k-1,M);a=[1:M];%%KINDA force the periodic options in chosiing the cell (Dleaware wave diffusion)

ii=[ii;p(a)]; jj=[jj;q(a)];%gain from the neigborh cell
end
ii=[ii;p];jj=[jj;p];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ACTUAL LOOP
for i=N-1:-1:1
        DD=(cg(i,:).*c(i,:))';
        %DD=c(i,:)';
        s=[];S=0*p;
        %%%%%%%%%%%%%%%%%%
        for k=[-1 1];
            
%         if periodic==0;     q=p+k;if k==1;a=[1:M-1];elseif k==-1;a=[2:M];end 
%         elseif periodic==1; q=1+mod(p+k-1,M);a=[1:M];
%         end
q=1+mod(p+k-1,M);a=[1:M];%%KINDA force the periodic options in chosiing the cell (Dleaware wave diffusion)
        
        %lateral diffraction
        %value=alpha/dx^2 /(2*(2*pi/T)) *((DD(p(a))+DD(q(a)))/2-1/2*DD(p(a))) *dx*stretchdx./cg(i+1,p(a))';   
        %value=alpha/dx^2 /(2*(2*pi/T)) *(DD(p(a))+DD(q(a)))/2 *dx*stretchdx./cg(i+1,p(a))';   
        value=alpha/dx^2 /(2*(2*pi/T)).*min(DD(p(a)),DD(q(a))) *dx*stretchdx./cg(i+1,p(a))'; %GOOD     
        %value=alpha/dx^2 /(2*(2*pi/T)) *min(DD(p(a)),DD(q(a))) *dx*stretchdx;   
        %value=alpha/dx^2 /(2*(2*pi/T)) *min(DD(p(a)),DD(q(a))) *dx*stretchdx./cg(i+1,p(a))';
        %value=alpha/dx^2 ./max(kwave(i,p(a)),kwave(i,q(a)))' *dx*stretchdx;   
        %value=alpha/dx^2 ./kwave(i+1,p(a))' *dx*stretchdx;  
           
            if periodic==0;%DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
            %REMOVE THE LATERAL DIFFUSION
            value=value.*(AW(i,p(a))==0)';%zero diffusion at the lateral boundary if not periodic
            end
        
        S(p(a))=S(p(a))+value; %exit from that cell
        s=[s;-value]; %gain from the neigborh cell
        end

        
        %Ei=(H(i+1,p).^2/16);  Ei is E!!!!
        %if abs(sum(E-Ei))>0.00001;pause;end
        
        %dissipation by breaking
        Qbrk=Qbforbreaking(8*(E./max(0.001,Hcr(i,p).^2)));
        brk=0.25/T*Qbrk.*Hcr(i,p).^2./E;
        brk(E<0.001)=0;
        
        
        %Whitecapping
        whitecap=0*3.33*10^-5*sigma.*(E*sigma^4/9.81^2/0.00457).^2;
        %Br=0.00175;
        %B=cg(i,p).*kwave(i,p).^3.*E;
        %whitecap=1000*5*10^-5*sqrt(9.81*kwave(i,p)).*(B/Br).^(2/2).*(B>Br);
        
        
        %Shoaling, bed friction, and breaking
        %s=[s;S+(  (cg(i,p)+dx*stretchdx*(bedfr(i,p))+brk+whitecap)./cg(i+1,p))'];
        s=[s;S+(  (cg(i,p)+dx*stretchdx*(bedfr(i,p)))./cg(i+1,p))'];


        A = sparse(ii,jj,s);%solve the matrix inversion
        E=(A\E')';
        
        %Set the lateral bouundaries if not periodics. Two options:
        %no-gradient (+-1) or zero wave heigth (+-2)
        if periodic==0 %DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
            if angle0>0%angle<0
                a=find(AW(i,:)==-1); %find the incomign boundary
                if a>0;
                b=find(AW(i,1+mod(-1+a-1,M))==0);
                E(a)=E(a(b)); 
                end   
                a=find(AW(i,:)==-2); %find the incomign boundary
                if a>0;
                b=find(AW(i,1+mod(-1+a-1,M))==0); %the first free cel;
                E(a)=E(a(b))*0+Holateral^2/16; 
                end 
            elseif angle0<0%angle>0
                a=find(AW(i,:)==1); %find the incomign boundary
                if a>0;
                b=find(AW(i,1+mod(-1+a+1,M))==0);
                E(a)=E(a(b)); 
                end   
                a=find(AW(i,:)==2); %find the incomign boundary
                if a>0;
                b=find(AW(i,1+mod(-1+a+1,M))==0);
                E(a)=E(a(b))*0+Holateral^2/16;  
                end  
            end
        end
    
        %recaluclate the wave height
        H(i,:)=4*sqrt(max(0,E));
        
        %breaking checks
        a=find(H(i,:)>Hcr(i,:));    
        H(i,a)=Hcr(i,a);
        E(a)=Hcr(i,a).^2/16;
        

        
end;%end of wave propagation

Uwave=(pi*H./T./sinh(kwave.*d));
Uwave(d<=0)=0;



%shear the domain
for i=1:N
    s=floor(tgangle*i);
    H(i,:)=circshift(H(i,:),[0 -s*angledir]);
    Uwave(i,:)=circshift(Uwave(i,:),[0 -s*angledir]);
end

%flip the domain
if gridDIR==-1; 
    H=H(end:-1:1,:);
    Uwave=Uwave(end:-1:1,:);
end

%figure;imagesc(H);pause
