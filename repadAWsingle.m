% function     AW=repadAWsingle(AW,angle0,N,M,i);
%     
%     if angle0>0 %wave right to left. Use the right boundary
%         a1=find(AWup(:)==-1 | AWup(:)==-2);
%         a2=find(AW(:)==-1 | AW(:)==-2);
%         if length(a2)>0  %this line introduced to put AW not along all side boundary, only a piece and allow some AW=0 along sides
%         if a2>=a1;AW(a1:a2)=AW(a2(1));else;AW(1:a2)=AW(a2(1));AW(a1:M)=AW(a2(1));end
%         end
%         
%     elseif angle0<0 %wave left to right. Use the left boundary
%         a2=find(AWup(:)==1 | AWup(:)==2);
%         a1=find(AW(:)==1 | AW(:)==2);
%         if length(a1)>0
%         if a2>=a1;AW(a1:a2)=AW(a1(1));else;AW(1:a2)=AW(a1(1));AW(a1:M)=AW(a1(1));end
%         end
%     end
%     

function     AW=repadAWsingle(AW,AWup,angle0,N,M,i);
    
    if angle0>0 %wave right to left. Use the right boundary
        a1=find(AWup(:)==-1 | AWup(:)==-2);
        a2=find(AW(:)==-1 | AW(:)==-2);
        if length(a2)>0  %this line introduced to put AW not along all side boundary, only a piece and allow some AW=0 along sides
        if a2>=a1;AW(a1:a2)=AW(a2(1));else;AW(1:a2)=AW(a2(1));AW(a1:M)=AW(a2(1));end
        end
        
    elseif angle0<0 %wave left to right. Use the left boundary
        a2=find(AWup(:)==1 | AWup(:)==2);
        a1=find(AW(:)==1 | AW(:)==2);
        if length(a1)>0
        if a2>=a1;AW(a1:a2)=AW(a1(1));else;AW(1:a2)=AW(a1(1));AW(a1:M)=AW(a1(1));end
        end
    end
    