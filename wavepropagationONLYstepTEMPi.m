function [E,dPW]=wavepropagationONLYstepTEMPi(gridDIR,E,Etot,EEE,i,cg,cgup,d,dup,bedfr,sigma,T,dx,angle,Cbr,Cbed,alpha,N,M,A,AW,AWup,p,ii,jj,periodic,Holateral,sigmatot)

PWup=E.*cgup;
%Shear the inputs
angledir=-sign(angle*gridDIR)*gridDIR;
tgangle=tan(abs(angle)/180*pi);

% s=floor(tgangle*(i+1));   
% E=circshift(E,[0 s*angledir]);   
% Etot=circshift(Etot,[0 s*angledir]);  
% EEE=circshift(EEE,[0 s*angledir]);
% cgup=circshift(cgup,[0 s*angledir]); 
% AWup=circshift(AWup,[0 s*angledir]); 
% dup=circshift(dup,[0 s*angledir]); 
% 
% s=floor(tgangle*(i));    
% sigma=circshift(sigma,[0 s*angledir]);
% sigmatot=circshift(sigmatot,[0 s*angledir]);  
% d=circshift(d,[0 s*angledir]);      
% bedfr=circshift(bedfr,[0 s*angledir]);      
% cg=circshift(cg,[0 s*angledir]);   
% AW=circshift(AW,[0 s*angledir]);   

%the cell before the step (upwind). the sign +1 is becuase we around
%counting i from N to -1 to 1 !!!. Thefore i+1 is more seaward than i!!!
s=floor(tgangle*(i+1));   
E=circshift(E,[0 s*angledir]);   
Etot=circshift(Etot,[0 s*angledir]);  
EEE=circshift(EEE,[0 s*angledir]);
cgup=circshift(cgup,[0 s*angledir]); 
AWup=circshift(AWup,[0 s*angledir]); 
dup=circshift(dup,[0 s*angledir]);d=dup;   

%the cell at the step
s=floor(tgangle*(i));       
cg=circshift(cg,[0 s*angledir]);   
AW=circshift(AW,[0 s*angledir]);   
%d=circshift(d,[0 s*angledir]);      
bedfr=circshift(bedfr,[0 s*angledir]);  
sigma=circshift(sigma,[0 s*angledir]);
sigmatot=circshift(sigmatot,[0 s*angledir]);     

AW=repadAWsingle(AW,AWup,angle,N,M,i); 


% angledir=-sign(angle*gridDIR)*gridDIR;
% tgangle=tan(abs(angle*gridDIR)/180*pi);
% 
% s=floor(tgangle*(i));
% AWup=circshift(AWup,[0 s*angledir]); 
% 
% s=floor(tgangle*(i+1));
% E=circshift(E,[0 s*angledir]);   
% Etot=circshift(Etot,[0 s*angledir]);  
% EEE=circshift(EEE,[0 s*angledir]);
% sigma=circshift(sigma,[0 s*angledir]);
% sigmatot=circshift(sigmatot,[0 s*angledir]);  
% cg=circshift(cg,[0 s*angledir]);   
% cgup=circshift(cgup,[0 s*angledir]);  
% AW=circshift(AW,[0 s*angledir]);   
% d=circshift(d,[0 s*angledir]);      
% bedfr=circshift(bedfr,[0 s*angledir]);   
%    
% 
% AW=repadAWsingle(AW,AWup,angle,N,M,i); 




%when you stretch the coordinates, you also stretch the dx!!!
stretchdx=1/cos(angle/180*pi);


        %Dissipation by breaking
        Hcr=max(0.01,Cbr*d(p));
        brk=0*E;
        val=min(1,(8*EEE./Hcr.^2));%EEE e' il totale, abdrebbe usato in teoria
        Qbrk=((val(val>0.1 & E>0.00001)-0.1)/0.9).^2;
        brk(val>0.1 & E>0.00001)=0.25/T*Qbrk.*(Hcr(val>0.1 & E>0.00001)).^2./E(val>0.1 & E>0.00001);%brk(E<0.01)=0;   

        %Whitecapping
        whitecap=3.33*10^-5*sigmatot.^7.*(EEE/9.81^2/0.00457.*sigma).^2;%
        
        %Bed friction
        bedfr=bedfr(p);
        
        
        %Shoaling, bed friction, and breaking
        E=E./((cg(p)+dx*stretchdx*(bedfr+whitecap+brk))./cgup(p));
        
        
        %Set the lateral bouundaries if not periodics. Two options:
        %no-gradient (+-1) or zero wave heigth (+-2)
        if periodic==0 %DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
            if angle>0%angle<0
                a=find(AW(:)==-1); %find the incomign boundary
                if a>0;
                b=find(AW(1+mod(-1+a-1,M))==0);
                E(a)=E(a(b)); 
                end   
                a=find(AW(:)==-2); %find the incomign boundary
                if a>0;
                b=find(AW(1+mod(-1+a-1,M))==0); %the first free cel;
                E(a)=E(a(b))*0+Holateral^2/16; 
                end 
            elseif angle<0%angle>0
                a=find(AW(:)==1); %find the incomign boundary
                if a>0;
                b=find(AW(1+mod(-1+a+1,M))==0);
                E(a)=E(a(b)); 
                end   
                a=find(AW(:)==2); %find the incomign boundary
                if a>0;
                b=find(AW(1+mod(-1+a+1,M))==0);
                E(a)=E(a(b))*0+Holateral^2/16;  
                end  
            end
        end
    
PW=E.*cg;
dPW=(PWup-PW)/(dx*stretchdx);

        
s=floor(tgangle*i); 
E=circshift(E,[0 -s*angledir]);
dPW=circshift(dPW,[0 -s*angledir]);
 