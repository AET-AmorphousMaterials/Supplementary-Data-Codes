function id=myunique(list)
list=sort(list,2);
% n=size(list,2);
% id=1;
% for i=2:size(list,1)
%     if ~sum(sum(list(id,:)==list(i,:),2)==n)
%         id=[id,i];
%     end
% end
[~, id] = unique(list, 'rows');
end