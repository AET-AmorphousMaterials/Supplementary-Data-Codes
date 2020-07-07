function [path3,path4,path5,path6,path7,path8] = uniqueCluster(path3,path4,path5,path6,path7,path8)
%% make sure the vertices are unique
path3(sum(path3==0,2)~=0,:)=[];
path4(sum(path4==0,2)~=0,:)=[];
path5(sum(path5==0,2)~=0,:)=[];
path6(sum(path6==0,2)~=0,:)=[];
path7(sum(path7==0,2)~=0,:)=[];
path8(sum(path8==0,2)~=0,:)=[];
if ~isempty(path3)
id1=sum((path3(:,2:3)-repmat(path3(:,1),[1 2]))==0,2)~=0;
id2=(path3(:,2)-path3(:,3))==0;
id=id1 | id2;
path3(id,:)=[];
end
if ~isempty(path4)
id1=sum((path4(:,2:end)-repmat(path4(:,1),[1 3]))==0,2)~=0;
id2=sum((path4(:,3:end)-repmat(path4(:,2),[1 2]))==0,2)~=0;
id3=(path4(:,4)-path4(:,3))==0;
id=id1 | id2 | id3;
path4(id,:)=[];
end
if ~isempty(path5)
id1=sum((path5(:,2:end)-repmat(path5(:,1),[1 4]))==0,2)~=0;
id2=sum((path5(:,3:end)-repmat(path5(:,2),[1 3]))==0,2)~=0;
id3=sum((path5(:,4:end)-repmat(path5(:,3),[1 2]))==0,2)~=0;
id4=(path5(:,5)-path5(:,4))==0;
id=id1 | id2 | id3 | id4;
path5(id,:)=[];
end
if ~isempty(path6)
id1=sum((path6(:,2:end)-repmat(path6(:,1),[1 5]))==0,2)~=0;
id2=sum((path6(:,3:end)-repmat(path6(:,2),[1 4]))==0,2)~=0;
id3=sum((path6(:,4:end)-repmat(path6(:,3),[1 3]))==0,2)~=0;
id4=sum((path6(:,5:end)-repmat(path6(:,4),[1 2]))==0,2)~=0;
id5=(path6(:,6)-path6(:,5))==0;
id=id1 | id2 | id3 | id4 | id5;
path6(id,:)=[];
end
if ~isempty(path7)
id1=sum((path7(:,2:end)-repmat(path7(:,1),[1 6]))==0,2)~=0;
id2=sum((path7(:,3:end)-repmat(path7(:,2),[1 5]))==0,2)~=0;
id3=sum((path7(:,4:end)-repmat(path7(:,3),[1 4]))==0,2)~=0;
id4=sum((path7(:,5:end)-repmat(path7(:,4),[1 3]))==0,2)~=0;
id5=sum((path7(:,6:end)-repmat(path7(:,5),[1 2]))==0,2)~=0;
id6=(path7(:,7)-path7(:,6))==0;
id=id1 | id2 | id3 | id4 | id5 | id6;
path7(id,:)=[];
end
if ~isempty(path8)
id1=sum((path8(:,2:end)-repmat(path8(:,1),[1 7]))==0,2)~=0;
id2=sum((path8(:,3:end)-repmat(path8(:,2),[1 6]))==0,2)~=0;
id3=sum((path8(:,4:end)-repmat(path8(:,3),[1 5]))==0,2)~=0;
id4=sum((path8(:,5:end)-repmat(path8(:,4),[1 4]))==0,2)~=0;
id5=sum((path8(:,6:end)-repmat(path8(:,5),[1 3]))==0,2)~=0;
id6=sum((path8(:,7:end)-repmat(path8(:,6),[1 2]))==0,2)~=0;
id7=(path8(:,8)-path8(:,7))==0;
id=id1 | id2 | id3 | id4 | id5 | id6 | id7;
path8(id,:)=[];
end
%% find unique polygons
if ~isempty(path3)
id=myunique(path3);
path3=path3(id,:);
end
if ~isempty(path4)
id=myunique(path4);
path4=path4(id,:);
end
if ~isempty(path5)
id=myunique(path5);
path5=path5(id,:);
end
if ~isempty(path6)
id=myunique(path6);
path6=path6(id,:);
end
if ~isempty(path7)
id=myunique(path7);
path7=path7(id,:);
end
if ~isempty(path8)
id=myunique(path8);
path8=path8(id,:);
end
%% remove path3 that share tetrahedron with path4-8
common=intersect(path3(:),path4(:));
if ~isempty(common)
id=sum(path3==common(1),2)~=0;
for i=common'
    id=id | sum(path3==i,2)~=0;
end
path3(id,:)=[];
end
common=intersect(path3(:),path5(:));
if ~isempty(common)
id=sum(path3==common(1),2)~=0;
for i=common'
    id=id | sum(path3==i,2)~=0;
end
path3(id,:)=[];
end
common=intersect(path3(:),path6(:));
if ~isempty(common)
id=sum(path3==common(1),2)~=0;
for i=common'
    id=id | sum(path3==i,2)~=0;
end
path3(id,:)=[];
end
common=intersect(path3(:),path7(:));
if ~isempty(common)
id=sum(path3==common(1),2)~=0;
for i=common'
    id=id | sum(path3==i,2)~=0;
end
path3(id,:)=[];
end
common=intersect(path3(:),path8(:));
if ~isempty(common)
id=sum(path3==common(1),2)~=0;
for i=common'
    id=id | sum(path3==i,2)~=0;
end
path3(id,:)=[];
end
end

