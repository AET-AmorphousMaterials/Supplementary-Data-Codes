function [path3,path4,path5,path6,path7,path8] = shareAB(path3,path4,path5,path6,path7,path8,connect)
%path3 will automatically share AB atom
% path4
ind=[];
for i=1:size(path4,1)
    list=connect(path4(i,:),:);
    a=intersect(list(1,:),list(2,:));
    b=intersect(list(3,:),list(4,:));
    if length(intersect(a,b))~=2
        ind=[ind i];
    end
end
path4(ind,:)=[];
% path5
ind=[];
for i=1:size(path5,1)
    list=connect(path5(i,:),:);
    a=intersect(list(1,:),list(2,:));
    b=intersect(list(3,:),list(4,:));
    b=intersect(b,list(5,:));
    if length(intersect(a,b))~=2
        ind=[ind i];
    end
end
path5(ind,:)=[];

% path6
ind=[];
for i=1:size(path6,1)
    list=connect(path6(i,:),:);
    a=intersect(list(1,:),list(2,:));
    b=intersect(list(3,:),list(4,:));
    c=intersect(list(5,:),list(6,:));
    b=intersect(b,c);
    if length(intersect(a,b))~=2
        ind=[ind i];
    end
end
path6(ind,:)=[];

% path7
ind=[];
for i=1:size(path7,1)
    list=connect(path7(i,:),:);
    a=intersect(list(1,:),list(2,:));
    a=intersect(a,list(7,:));
    b=intersect(list(3,:),list(4,:));
    c=intersect(list(5,:),list(6,:));
    b=intersect(b,c);
    if length(intersect(a,b))~=2
        ind=[ind i];
    end
end
path7(ind,:)=[];
% path8
ind=[];
for i=1:size(path8,1)
    list=connect(path8(i,:),:);
    a=intersect(list(1,:),list(2,:));
    b=intersect(list(3,:),list(4,:));
    c=intersect(list(5,:),list(6,:));
    d=intersect(list(7,:),list(8,:));
    a=intersect(a,b);
    b=intersect(c,d);
    if length(intersect(a,b))~=2
        ind=[ind i];
    end
end
path8(ind,:)=[];

end

