%% Main_3_polytetrahedral_analysis
% This script searches for polytetrahedral motifs
clear
clc
addpath('./src')

% perform analysis using different distortion parameters
for delta=[0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.255 0.275 0.3 0.325 0.35]
% read in files: finalized atomic coordinates in Angstrom and normalized
% BOO
data=importdata('./Input/Ta_film/BOO_normalized.mat');
model=data.model;
order=data.order;

% separate crystalline and amorphus atom by using BOO at 0.5
ind1=order>0.5;
ind2=~ind1;

% remove scattered crystalline atoms that are isolated from others
atom=ones(1,size(model,2)); atom(ind1)=1;   atom(ind2)=2;
model=double(model);

ind1=find(ind1);
ind2=find(ind2);
temp=[];
for i=1:length(ind1)
    dis=model(:,ind1)-model(:,ind1(i));
    dis=sqrt(sum(dis.^2,1));
    indt=find(dis<=8);
    if length(indt)<13
        temp=[temp i];
    end
end
ind3=ind1(temp);    ind1(temp)=[];

temp=[];
for i=1:length(ind3)
    dis=model(:,ind1)-model(:,ind3(i));
    dis=sqrt(sum(dis.^2,1));
    indt=find(dis<=8);
    if length(indt)>1
        temp=[temp i];
    end
end
ind1=[ind1 ind3(temp)]; ind3(temp)=[];
ind2=[ind2 ind3];

modelcrystal=model(:,ind1);
modelamorphous=model(:,ind2);
model=double([modelamorphous,modelcrystal]);
DT = delaunayTriangulation(model(1,:)',model(2,:)',model(3,:)');

% remove tetrahedron with center outside the sample
DTremove=DT;
IC = incenter(DTremove);
shp = alphaShape(model(1,:)',model(2,:)',model(3,:)',4);

in = inShape(shp,IC(:,1),IC(:,2),IC(:,3));
DTremove = triangulation(DTremove.ConnectivityList(in,:),model');

% Identify acceptable tetrahedra by using distortion parameter and 1st
% nearest shell distance
ConnectivityList=DTremove.ConnectivityList;
ind=[];
for i=1:size(ConnectivityList,1)
    points=model(:,ConnectivityList(i,:));
    edge=[points(:,2:4)-points(:,1), points(:,3:4)-points(:,2), points(:,4)-points(:,3)];
    edge=sqrt(sum(edge.^2,1));
    if max(edge)>mean(edge)*(delta+1) | max(edge)>4.12 % 1st valley in RDF
        ind=[ind,i];
    end
end
ConnectivityList(ind,:)=[];
DTremove = triangulation(ConnectivityList,model');
IC = incenter(DTremove);

% find face sharing tetrahedra
Ntet=size(DTremove.ConnectivityList,1);
Neiface=[];
for i=1:Ntet-1
    i
    vertex1=DTremove.ConnectivityList(i,:);
    temp=DTremove.ConnectivityList==vertex1(1) | DTremove.ConnectivityList==vertex1(2) | DTremove.ConnectivityList==vertex1(3) | DTremove.ConnectivityList==vertex1(4);
    temp=sum(temp,2);
    jrange=find(temp~=0);
    jrange(jrange<=i)=[];
    if ~isempty(jrange)
    for j=jrange'
        %fast methd to find number of common elements
        com = DTremove.ConnectivityList([i j],:);
        com=sort(com(:));
        temp=sum(com([2 4 6])==com([1 3 5]))+sum(com([2 4 6])==com([3 5 7]))+sum(com(7)==com(8));
        if temp==3 % share 3 vertex
            Neiface=[Neiface;i j];
        end
    end
    end
end
% find different 3 to 8 fold motifs
[path3,path4,path5,path6,path7,path8]=findclustersN(IC,Neiface);
% find unique motifs and remove 3-fold motifs that share atoms with
% other types of motifs
[path3,path4,path5,path6,path7,path8] = ...
    uniqueCluster(path3,path4,path5,path6,path7,path8);
% remove large motifs that contain smaller motifs
[path3,path4,path5,path6,path7,path8] = smallestCluster(path3,path4,path5,path6,path7,path8,Neiface);
% all motifs must share same AB atom
[path3,path4,path5,path6,path7,path8] = shareAB(path3,path4,path5,path6,path7,path8,DTremove.ConnectivityList);
% save results
eval(['save ./Output/Ta_film/Polytetrahedral_motifs_delta_' num2str(delta) '.mat'])
end

% plot results
results=[]; ni=0;
for delta=[0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.255 0.275 0.3 0.325 0.35]
ni=ni+1;
eval(['load ./Output/Ta_film/Polytetrahedral_motifs_delta_' num2str(delta) '.mat']) % laod above analysis results
amorInd=1:size(modelamorphous,2);  % results for amorphous atoms
ICamor=find(sum(DTremove.ConnectivityList<=amorInd(end) & DTremove.ConnectivityList>=amorInd(1),2)>2);

temp=intersect(unique(DTremove.ConnectivityList(:)),amorInd);
results(1,ni)=length(temp)/length(amorInd);

% fraction based on tetrahedra
results(2,ni)=length(unique(intersect(path3(:),ICamor)))/size(ICamor,1);
results(3,ni)=length(unique(intersect(path4(:),ICamor)))/size(ICamor,1);
results(4,ni)=length(unique(intersect(path5(:),ICamor)))/size(ICamor,1);
results(5,ni)=length(unique(intersect(path6(:),ICamor)))/size(ICamor,1);
results(6,ni)=length(unique(intersect(path7(:),ICamor)))/size(ICamor,1);
results(7,ni)=length(unique(intersect(path8(:),ICamor)))/size(ICamor,1);
end
% plot polytetrahedral analysis results as a function of delta
figure(13)
h = plot(1:10,zeros(10,10));
c = get(h,'Color'); % get colors
temp=c{1};
c{1}=c{4};
c{4}=c{3};
c{3}=c{2};
c{2}=temp;
beta=[0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.255 0.275 0.3 0.325 0.35]';
Ta=results'*100;

figure(13)
clf
plot(beta,Ta(:,2),'o-','linewidth',2,'Color',c{1})      
hold on
plot(beta,Ta(:,3),'o-','linewidth',2,'Color',c{2}) 
plot(beta,Ta(:,4),'o-','linewidth',2,'Color',c{3}) 
plot(beta,Ta(:,5),'o-','linewidth',2,'Color',c{4}) 
% plot(beta,Ta(:,6),'o-','linewidth',2,'Color',c{5}) 
% plot(beta,Ta(:,7),'o-','linewidth',2,'Color',c{6}) 
plot(beta,Ta(:,1),'o-','linewidth',2,'Color',c{5}) 
xlim([0 0.35])
set(gca,'xTick',0:0.05:0.35,'xticklabel',0:0.05:0.35)
hold off
box off
 ax = axis;
 hold on
 plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1)
 plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1)
 hold off
ylabel('Fraction (%)')
xlabel('\delta')
set(gca,'fontsize', 16,'FontName', 'Arial')
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth=1;
set(gca,'FontSize',16,'FontName', 'Arial','fontweight','bold');
legend('Triplet', 'Quadrilateral','Pentagon','Hexagon','Tetrahedral atom','FontSize',14,'fontweight','normal')%
legend boxoff