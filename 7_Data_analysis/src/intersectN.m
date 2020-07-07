function [ele] = intersectN(data)
ele=intersect(data(1,:),data(2,:));
if size(data,1)>2
for i=3:size(data,1)
    ele=intersect(ele,data(i,:));
end
end
end

