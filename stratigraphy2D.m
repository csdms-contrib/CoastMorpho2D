function [Y1,Y2,Yb,flyr,flyrb,plyr]=stratigraphy2D(A,Y1,Y2,Yb,flyr,flyrb,plyr,nlyr,dlyr,tlyrU,tlyrD);

%to temove the negatives
StY1=max(0,-Y1);StY2=max(0,-Y2);
Y1=Y1+StY1;Y2=Y2+StY2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Y1+Y2;

exsE=find(Y>tlyrU & A==1); %accretion past the limit
exsD=find(Y<tlyrD & A==1); %erosion past the limit
    


%accretion past limit -> need to remove material from top cell
    flyrU=Y1./(Y1+Y2);    
    deltaY1=dlyr*flyrU(exsE);%how much is removed
    deltaY2=dlyr*(1-flyrU(exsE));%how much is removed
    Y1(exsE)=Y1(exsE)-deltaY1;%remove it from the active layer
    Y2(exsE)=Y2(exsE)-deltaY2;%remove it from the active layer
   
    
    %add the stuff back to the top layer         
        ao=find(plyr(exsE)<nlyr);
        a=find(plyr(exsE)==nlyr);
        
        %add a new layer on top
        plyr(exsE(ao))=plyr(exsE(ao))+1;
        
        %shift the layer compostion, so that you can fill in the top layer
        Ybnew=Yb(exsE(a))+dlyr;%the new thickness of the bottom layer
        firstlyr=squeeze(flyr(:,:,1));flyrb(exsE(a))=(flyrb(exsE(a)).*Yb(exsE(a))+firstlyr(exsE(a))*dlyr)./Ybnew;%store the fraction in the bottom layer
        Yb(exsE(a))=Ybnew;%increase the thickness of the bed layer
        for i=1:nlyr-1;
        temp=squeeze(flyr(:,:,i));temp2=squeeze(flyr(:,:,i+1));
        temp(exsE(a))=temp2(exsE(a));
        flyr(:,:,i)=temp;%make the actual shift in layer composition
        end
        
        %eiether way, update the top layer  
        for i=1:nlyr
            temp=squeeze(flyr(:,:,i));
            a=find(plyr(exsE)==i);temp(exsE(a))=flyrU(exsE(a));
            flyr(:,:,i)=temp;
        end
        
        
        
        %YOU CAN ONLY TAKE SEDIEMNT FROM ONE LAYER AT THE TIME
%erosion past the limit: need to take from bottom cell

        a=find(plyr(exsD)>=1); %the level occupied is above the base layer
        %plyr(exsD)
        %pause
        ao=find(plyr(exsD)==0);
        
        deltaY1=0*Y;deltaY2=0*Y;
        
        %take it from the uppermost layer
        for i=1:nlyr;            %go over all the slices
            temp=squeeze(flyr(:,:,i));
            aa=find(plyr(exsD(a))==i);
            deltaY1(exsD(a(aa)))=temp(exsD(a(aa)))*dlyr;
            deltaY2(exsD(a(aa)))=(1-temp(exsD(a(aa))))*dlyr;
            
            temp(exsD(a(aa)))=NaN;%clear the cell that you eliminated
            flyr(:,:,i)=temp;%clear the cell that you eliminated
        end
        plyr(exsD(a))=plyr(exsD(a))-1;%lower the elevation of the top layer
        
       
        %take it from the bottom layer (if your at the lowest level)
        deltaYb=min(dlyr,Yb(exsD(ao)));
        deltaY1(exsD(ao))=deltaYb.*flyrb(exsD(ao));
        deltaY2(exsD(ao))=deltaYb.*(1-flyrb(exsD(ao)));
        Yb(exsD(ao))=Yb(exsD(ao))-deltaYb;
        
        
        %either way (normal layer or bottom layer), add it!, becuase you already removed it
        Y1=Y1+deltaY1;%how much is added
        Y2=Y2+deltaY2;%how much is added
        
       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%re-add the negative layers
Y1=Y1-StY1;Y2=Y2-StY2;