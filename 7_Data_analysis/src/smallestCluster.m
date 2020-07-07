function [path3,path4,path5,path6,path7,path8] = smallestCluster(path3,path4,path5,path6,path7,path8,Neiface)
%% remove large polygons that contain smaller polygons
temp=sort(Neiface,2);
temp=temp(:,1)+temp(:,2)*1e7;
% squres
if ~isempty(path4)
temp1=sort(path4(:,[1 3]),2);
temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path4(:,[2 4]),2);
temp2=temp2(:,1)+temp2(:,2)*1e7;
id=zeros(size(path4,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
if sum(id)~=0
path4(id,:)=[];
end
end
% pentagon
if ~isempty(path5)
temp1=sort(path5(:,[1 3]),2);
temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path5(:,[2 4]),2);
temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path5(:,[3 5]),2);
temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path5(:,[4 1]),2);
temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path5(:,[5 2]),2);
temp5=temp5(:,1)+temp5(:,2)*1e7;
id=zeros(size(path5,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
if sum(id)~=0
path5(id,:)=[];
end

% hexagon
if ~isempty(path6)
temp1=sort(path6(:,[1 4]),2);
temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path6(:,[2 5]),2);
temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path6(:,[3 6]),2);
temp3=temp3(:,1)+temp3(:,2)*1e7;
id=zeros(size(path6,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
if sum(id)~=0
path6(id,:)=[];
end

temp1=sort(path6(:,[1 3]),2);
temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path6(:,[2 4]),2);
temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path6(:,[3 5]),2);
temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path6(:,[4 6]),2);
temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path6(:,[5 1]),2);
temp5=temp5(:,1)+temp5(:,2)*1e7;
temp6=sort(path6(:,[6 2]),2);
temp6=temp6(:,1)+temp6(:,2)*1e7;
id=zeros(size(path6,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
c=intersect(temp6,temp);
for i=c'
    id=id | temp6==i;
end
if sum(id)~=0
path6(id,:)=[];
end
end

% octagon
if ~isempty(path8)
temp1=sort(path8(:,[1 4]),2);temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path8(:,[2 5]),2);temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path8(:,[3 6]),2);temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path8(:,[4 7]),2);temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path8(:,[5 8]),2);temp5=temp5(:,1)+temp5(:,2)*1e7;
temp6=sort(path8(:,[6 1]),2);temp6=temp6(:,1)+temp6(:,2)*1e7;
temp7=sort(path8(:,[7 2]),2);temp7=temp7(:,1)+temp7(:,2)*1e7;
temp8=sort(path8(:,[8 3]),2);temp8=temp8(:,1)+temp8(:,2)*1e7;
id=zeros(size(path8,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
c=intersect(temp6,temp);
for i=c'
    id=id | temp6==i;
end
c=intersect(temp7,temp);
for i=c'
    id=id | temp7==i;
end
c=intersect(temp8,temp);
for i=c'
    id=id | temp8==i;
end
if sum(id)~=0
path8(id,:)=[];
end

temp1=sort(path8(:,[1 3]),2);temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path8(:,[2 4]),2);temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path8(:,[3 5]),2);temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path8(:,[4 6]),2);temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path8(:,[5 7]),2);temp5=temp5(:,1)+temp5(:,2)*1e7;
temp6=sort(path8(:,[6 8]),2);temp6=temp6(:,1)+temp6(:,2)*1e7;
temp7=sort(path8(:,[7 1]),2);temp7=temp7(:,1)+temp7(:,2)*1e7;
temp8=sort(path8(:,[8 2]),2);temp8=temp8(:,1)+temp8(:,2)*1e7;
id=zeros(size(path8,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
c=intersect(temp6,temp);
for i=c'
    id=id | temp6==i;
end
c=intersect(temp7,temp);
for i=c'
    id=id | temp7==i;
end
c=intersect(temp8,temp);
for i=c'
    id=id | temp8==i;
end
if sum(id)~=0
path8(id,:)=[];
end

temp1=sort(path8(:,[1 5]),2);temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path8(:,[2 6]),2);temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path8(:,[3 7]),2);temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path8(:,[4 8]),2);temp4=temp4(:,1)+temp4(:,2)*1e7;
id=zeros(size(path8,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
if sum(id)~=0
path8(id,:)=[];
end
end

% septagon
if ~isempty(path7)
temp1=sort(path7(:,[1 4]),2);temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path7(:,[2 5]),2);temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path7(:,[3 6]),2);temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path7(:,[4 7]),2);temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path7(:,[5 1]),2);temp5=temp5(:,1)+temp5(:,2)*1e7;
temp6=sort(path7(:,[6 2]),2);temp6=temp6(:,1)+temp6(:,2)*1e7;
temp7=sort(path7(:,[7 3]),2);temp7=temp7(:,1)+temp7(:,2)*1e7;
id=zeros(size(path7,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
c=intersect(temp6,temp);
for i=c'
    id=id | temp6==i;
end
c=intersect(temp7,temp);
for i=c'
    id=id | temp7==i;
end
if sum(id)~=0
path7(id,:)=[];
end

temp1=sort(path7(:,[1 3]),2);temp1=temp1(:,1)+temp1(:,2)*1e7;
temp2=sort(path7(:,[2 4]),2);temp2=temp2(:,1)+temp2(:,2)*1e7;
temp3=sort(path7(:,[3 5]),2);temp3=temp3(:,1)+temp3(:,2)*1e7;
temp4=sort(path7(:,[4 6]),2);temp4=temp4(:,1)+temp4(:,2)*1e7;
temp5=sort(path7(:,[5 7]),2);temp5=temp5(:,1)+temp5(:,2)*1e7;
temp6=sort(path7(:,[6 1]),2);temp6=temp6(:,1)+temp6(:,2)*1e7;
temp7=sort(path7(:,[7 2]),2);temp7=temp7(:,1)+temp7(:,2)*1e7;
id=zeros(size(path7,1),1);
c=intersect(temp1,temp);
for i=c'
    id=id | temp1==i;
end
c=intersect(temp2,temp);
for i=c'
    id=id | temp2==i;
end
c=intersect(temp3,temp);
for i=c'
    id=id | temp3==i;
end
c=intersect(temp4,temp);
for i=c'
    id=id | temp4==i;
end
c=intersect(temp5,temp);
for i=c'
    id=id | temp5==i;
end
c=intersect(temp6,temp);
for i=c'
    id=id | temp6==i;
end
c=intersect(temp7,temp);
for i=c'
    id=id | temp7==i;
end
if sum(id)~=0
path7(id,:)=[];
end
end
end

