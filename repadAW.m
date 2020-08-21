function     AW=repadAW(AW,angle0,N,M);
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