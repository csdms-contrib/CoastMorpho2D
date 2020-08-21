function [E]=wavepropagationstep(E,i,c,cg,d,kwave,sigma,T,dx,angle,Cbr,Cbed,alpha,angle0,N,M,A,AW,p,ii,jj,periodic,Holateral)

%when you stretch the coordinates, you also stretch the dx!!!
stretchdx=1/cos(angle/180*pi);
%wave breaking crietion
%Hcr=0.14*2*pi./kwave.*tanh(kwave.*d);Hcr(d<=0.1)=0;
%Hcr=Cbr*d;Hcr(d<=0.1)=0.01;



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
        Hcr=max(0.01,Cbr*d(i,p));
        Qbrk=Qbforbreaking(8*(E./max(0.001,Hcr.^2)));
        brk=0.25/T*Qbrk.*Hcr.^2./E;
        brk(E<0.001)=0;
        
        
        %Whitecapping
        whitecap=3.33*10^-5*sigma.*(E*sigma^4/9.81^2/0.00457).^2;
        %Br=0.00175;
        %B=cg(i,p).*kwave(i,p).^3.*E;
        %whitecap=1000*5*10^-5*sqrt(9.81*kwave(i,p)).*(B/Br).^(2/2).*(B>Br);
        
        %bed friction
        bedfr=Cbed*(sigma./(9.81*sinh(kwave(i,p).*d(i,p)))).^2;bedfr(d(i,p)<=0.1)=1;
        %bedfr=0.067*(sigma./(9.81*sinh(kwave.*d))).^2;bedfr(d<=0.1)=1;
        
        
        %Shoaling, bed friction, and breaking
        s=[s;S+(  (cg(i,p)+dx*stretchdx*(bedfr+brk+whitecap))./cg(i+1,p))'];
        %s=[s;S+(  (cg(i,p)+dx*stretchdx*(bedfr(i,p))+brk+whitecap)./cg(i+1,p))'];


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
    
 