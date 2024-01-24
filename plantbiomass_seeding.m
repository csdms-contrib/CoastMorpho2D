function [bio]=plantbiomass_seeding(A,dx,dt,bioseed,z,bio,Zlev,dBseed);


p=find(A==1); %how many points can become a pond, how much normal cells are present.

a=find(bio==0 & Zlev>dBseed);%only create new ponds when the scour can be made!

r=rand(length(a),1);

aa=a(r<bioseed*dt/365);
for i=1:length(aa)
    bio(aa(i))=0.1;
end

