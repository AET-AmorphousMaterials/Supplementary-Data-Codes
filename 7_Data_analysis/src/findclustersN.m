function [path3,path4,path5,path6,path7,path8] = findclustersN(IC,Neiface)
path3=cell(size(IC,1),1);
path4=cell(size(IC,1),1);
path5=cell(size(IC,1),1);
path6=cell(size(IC,1),1);
path7=cell(size(IC,1),1);
path8=cell(size(IC,1),1);

parfor node1=1:size(IC,1)
    node1
    chain=node1;
    for i=2:5 %search for connected tetrahera
        chaintemp=[];
        for nj=1:length(chain(:,i-1))
            j=chain(nj,i-1);
            temp=sum(Neiface==j,2)~=0;
            temp=Neiface(temp,:);
            temp=unique(temp(:));
            temp(temp==j)=[];
            if i>2
                temp(temp==chain(nj,i-2))=[];
            end
            if isempty(temp)
                temp=0; %if no next neighbor found, set the next one to be zero
            end
            chaintemp=[chaintemp; ones(length(temp),1)*chain(nj,:) temp];
        end
        chain=chaintemp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find broken pentagon
    % both end are broken
%     id=chain(:,3)==0;
%     if sum(id)>=2
%     common1=unique(chain(id,2));
%     path=[];
%     for i=1:length(common1)-1
%         for j=i+1:length(common1)
%             path=[path;common1(i) chain(1,1) common1(j)];
%         end
%     end
%     path3=[path3;path];
%     end

    % at least one end is broken
    id=chain(:,3)~=0 & chain(:,4)==0;
    path3{node1}=chain(id,1:3);
%     path3=[path3; chain(id,1:3)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find squares
    common1=unique(chain(:,3)); 
    [yy,common1]=hist(chain(:,3),common1);
    common1=common1(yy>1);
    if ~isempty(common1)
        path4t=[];
    for common=common1'
        if common==0
            continue
        end
    id=chain(:,3)==common;
    path1=chain(id,1:3);
    id=myunique(path1);
    path1=path1(id,:);
    path=[];
    for i=1:size(path1,1)-1
        for j=i+1:size(path1,1)
            path=[path;path1(i,:) path1(j,2)];
        end
    end
    path4t=[path4t;path];    
    end
    path4{node1}=path4t;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common1=intersect(chain(:,3),chain(:,4)); %find pentagon
    if ~isempty(common1)
        path5t=[];
    for common=common1'
        if common==0
            continue
        end
    id=chain(:,4)==common;
    path1=chain(id,1:4);
    id=myunique(path1);
    path1=path1(id,:);
    
    id=chain(:,3)==common;
    path2=chain(id,2);
    path2=unique(path2); % chains shorter than the maxium length needs to be unique (the maxium chain is auto guaranteed)
    path=[];
    for i=1:length(path2)
        path=[path;path1 ones(size(path1,1),1)*path2(i)];
    end
    path5t=[path5t;path];
    end
    path5{node1}=path5t;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common1=unique(chain(:,4)); %find hexagon
    [yy,common1]=hist(chain(:,4),common1);
    common1=common1(yy>1);
    if ~isempty(common1)
        path6t=[];
    for common=common1'
        if common==0
            continue
        end
    id=chain(:,4)==common;
    path1=chain(id,1:4);
    id=myunique(path1);
    path1=path1(id,:);
    path=[];
    for i=1:size(path1,1)-1
        for j=i+1:size(path1,1)
            path=[path;path1(i,:) path1(j,3:-1:2)];
        end
    end
    path6t=[path6t;path];    
    end
    path6{node1}=path6t;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common1=intersect(chain(:,4),chain(:,5)); %find heptagon
    if ~isempty(common1)
        path7t=[];
    for common=common1'
        if common==0
            continue
        end
    id=chain(:,5)==common;
    path1=chain(id,1:5);
    id=myunique(path1);
    path1=path1(id,:);
    
    id=chain(:,4)==common;
    path2=chain(id,2:3);
    id=myunique(path2);
    path2=path2(id,:);
    path=[];
    for i=1:size(path2,1)
        path=[path;path1 ones(size(path1,1),1)*path2(i,2:-1:1)];
    end
    path7t=[path7t;path];
    end
    path7{node1}=path7t;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    common1=unique(chain(:,5)); %find octagon
    [yy,common1]=hist(chain(:,5),common1);
    common1=common1(yy>1);
    if ~isempty(common1)
        path8t=[];
    for common=common1'
        if common==0
            continue
        end
    id=chain(:,5)==common;
    path1=chain(id,1:5);
    id=myunique(path1);
    path1=path1(id,:);
    path=[];
    for i=1:size(path1,1)-1
        for j=i+1:size(path1,1)
            path=[path;path1(i,:) path1(j,4:-1:2)];
        end
    end
    path8t=[path8t;path];
    end
    path8{node1}=path8t;
    end
    
end
path3=cell2mat(path3);
path4=cell2mat(path4);
path5=cell2mat(path5);
path6=cell2mat(path6);
path7=cell2mat(path7);
path8=cell2mat(path8);
end

