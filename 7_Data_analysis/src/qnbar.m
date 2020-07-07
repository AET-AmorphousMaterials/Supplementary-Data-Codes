function [q] = qnbar(n,model,cutoff)
num=size(model,2);
q=zeros([1 num]);
qlm=zeros([2*n+1 num]);
% Qlm=zeros([2*n+1 num]);
for i=1:num
    dis=model(:,i)-model;
    dis=sqrt(dis(1,:).^2+dis(2,:).^2+dis(3,:).^2);
    ind=find(dis<cutoff & dis>0.1);
    if ~isempty(ind)
        vecs=model(:,ind)-model(:,i);
        [ph,th,~] = cart2sph(vecs(1,:),vecs(2,:),vecs(3,:));
        th=pi/2-th;
        for m=-n:n
            Y = harmonicY(n,m,th,ph);
            qlm(n+1+m,i)=sum(Y)/length(ind);
        end
    end
end
for i=1:num
    dis=model(:,i)-model;
    dis=sqrt(dis(1,:).^2+dis(2,:).^2+dis(3,:).^2);
    ind=find(dis<cutoff);% & dis>0.1);
    q(i)=sum(abs(sum(qlm(:,ind),2)/length(ind)).^2);
%     Qlm(:,i)=sum(qlm(:,ind),2)/max(1,length(ind));
end
q=sqrt(q*4*pi/(2*n+1));
% q=Qlm;
end

