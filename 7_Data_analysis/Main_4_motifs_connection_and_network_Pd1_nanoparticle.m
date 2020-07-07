%% Main_4_motifs_connection_and_network
% This script analyzes the edge-/vertex-sharing between 5-fold motifs,
% icosahedral filling by 5-fold motifs, as well as searches for pentagonal bipyramid
% networks
clear
clc
addpath('./src')

% load polytetrahedral analysis results with delta=0.255
load ./Input/Pd1_nanoparticle/Polytetrahedral_motifs_delta_0.255.mat

% select the motifs in the amorphous region
amorInd=1:size(modelamorphous,2);
ICamor=find(sum(DTremove.ConnectivityList<=amorInd(end) & DTremove.ConnectivityList>=amorInd(1),2)>2);
ind=[];
for i=1:size(path5,1)
    if intersect(path5(i,:),ICamor)>=length(path5(i,:))/2
        ind=[ind i];
    end
end
path5=path5(ind,:);

ind=[];
for i=1:size(path6,1)
    if intersect(path6(i,:),ICamor)>=length(path6(i,:))/2
        ind=[ind i];
    end
end
path6=path6(ind,:);

ind=[];
for i=1:size(path4,1)
    if intersect(path4(i,:),ICamor)>=length(path4(i,:))/2
        ind=[ind i];
    end
end
path4=path4(ind,:);

% search for edge-/vertex-sharing pentagonal bipyramids
Cor=[];Edg=[];
Cor=cell(15,1);
Edg=cell(15,1);
for i=1:size(path5,1)
    id=path5(i,:);
    temp=sum(path5==id(1) | path5==id(2) | path5==id(3) | path5==id(4) | path5==id(5),2);
    temp(i)=0;
    id1=find(temp~=0);
    corN=0;
    corind=[];
    edgN=0;
    edgind=[];
    for j=id1'
        temp=intersect(id,path5(j,:));
        if length(temp)==1
            corN=corN+1;
            corind=[corind j];
        elseif length(temp)==2
            edgN=edgN+1;
            edgind=[edgind j];
        end
    end
    Cor{corN+1}=[Cor{corN+1}; i corind];
    Edg{edgN+1}=[Edg{edgN+1}; i edgind];
end

% analyze icosahedron site filling
AB=[];
for i=1:size(path5,1)
    id=path5(i,:);
    AB(i,:)=intersectN(DTremove.ConnectivityList(id,:));
end
icosfil=[];
ABind=unique(AB(:));
icosfil=cell(12,1);
for i=ABind'
    temp=sum(AB==i,2);
    N=sum(temp);
    icosfil{N}=[icosfil{N}; find(temp==1)'];
end
icosfil{1}=unique(icosfil{1});% isolated pentagonal bipyramids can only contribute only 1 icosahedral site

% plot edge-/vertex-sharing pentagonal bipyramids results
tot=size(path5,1);
NN=[];
for i=1:size(Cor,1)
NN(1,i)=size(Cor{i},1);
end
for i=1:size(Edg,1)
NN(2,i)=size(Edg{i},1);
end
figure(1)
neifilling=NN'/tot;
neifilling(neifilling~=0)=neifilling(neifilling~=0)+0.003;
b=bar(neifilling*100,'BarWidth', 1);
b(1).FaceColor=[0, 0.4470, 0.7410];
b(2).FaceColor=[0.6350, 0.0780, 0.1840];
ylabel('Fraction (%)')
xlabel('Nerighboring pentagons')
set(gca,'xticklabel',0:14)
legend('Vertex sharing', 'Edge sharing')
legend boxoff
xlim([0.5 11.5])
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth=1.5;
set(gca,'FontSize',16,'FontName', 'Arial','fontweight','bold');

% plot icosahedron site filling
NN=[];
for i=1:size(icosfil,1)
NN(i)=size(icosfil{i},1);
end
tot=sum(NN);
figure(2)
icofil=NN'/tot;
icofil(icofil~=0)=icofil(icofil~=0)+0.003;
bar(icofil*100)
ylabel('Fraction (%)')
xlabel('Icosahedral filling')
set(gca,'xticklabel',1:12)
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth=1;
set(gca,'FontSize',16,'FontName', 'Arial','fontweight','bold');

%% Search for networks formed by quadrilateral/pentagonal/hexagonal motifs
% loop through 4/5/6-fold motifs
for iii=[4 5 6]
eval(['path=path' num2str(iii) ';']); % choose 4/5/6-fold motifs for network analysis

% find all edge sharing pairs
pair=[];
for i=1:size(path,1)-1
    i
    for j=i+1:size(path,1)
        if length(intersect(path(i,:),path(j,:)))==2
            pair=[pair;
                i j];
        end
    end
end

% search for networks of given motifs
connect=zeros(500,500);
nn=0;
for i=1:size(path,1)
    i
    temp=unique(connect(:));
    if sum(i==temp)==0
    list=cell(50,1);
    list{1}=i;
    n=1;
    while ~isempty(list{n})
        temp=pair==list{n}(1);
        for j=2:size(list{n},1)
            temp=temp | pair==list{n}(j);
        end
    ind=find(sum(temp,2)~=0);
    temp=pair(ind,:);
    temp=unique(temp(:));
    temp1=[];
        for j=1:n
            temp1=[temp1;list{j}];
        end
    list{n+1}=setdiff(temp,temp1);
    n=n+1;
    end
    temp1=[];
        for j=1:n
            temp1=[temp1;list{j}];
        end
    connect(nn+1,1:length(temp1))=temp1;
    nn=nn+1;
    end
end
save(['./Output/Pd1_nanoparticle/Network_of_' num2str(iii) '_fold_motifs.mat'])
end
% plot network size statistics
xx=5:1:200;
yy=[];
xxA=0:30;
yyA=[];
for i=4:6
    data=importdata(['./Output/Pd1_nanoparticle/Network_of_' num2str(i) '_fold_motifs.mat']);
    connectN=sum(data.connect~=0,2);
    connectN(connectN<5)=[];
    [y]=hist(connectN,xx);
    yy(i-3,:)=y;
    connect=data.connect;
    path=data.path;
    ConnectivityList=data.ConnectivityList;
    model=data.model;
    MROleg=[];
    for j=1:size(connect,1)
        if sum(connect(j,:)~=0)<5
            continue
        end
        ind=connect(j,connect(j,:)~=0);
        ind=path(ind,:);
        ind=unique(ind(:));
        ind=ConnectivityList(ind,:);
        ind=unique(ind(:));
        legtemp=0;
        for k=1:length(ind)-1
            for l=k+1:length(ind)
                dif=model(:,ind(k))-model(:,ind(l));
                dis=sqrt(sum(dif.^2));
                if dis>legtemp
                    legtemp=dis;
                end
            end
        end
        MROleg=[MROleg;legtemp];
    end
    [y]=hist(MROleg,xxA);
    yyA(i-3,:)=y;
end
s1=40;
ind=135;
ind1=find(xx==s1);
ind2=find(xx==ind);
xx(ind1+1:end)=[];
xx(ind1+1)=s1+1;
yy(:,ind1+1)=yy(:,ind2);
yy(:,ind1+2:end)=[];

% network size
figure(3)
clf
b=bar(xx,yy(:,:)','stacked');
b(1).FaceColor=[0, 0.4470, 0.7410];
ylim([0 30]);hold off;ylabel('Number of networks');xlabel('Network size')
set(gca,'xtick',[5:5:s1-5 s1+1],'xticklabel',[5:5:s1-5 ind])
box off;ax = axis;hold on
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1.5)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1.5);hold off
legend('Quadrilateral','Pentagonal','Hexagonal');legend boxoff
ax = gca;ax.LineWidth=1.5;set(gca,'FontSize',16,'FontName', 'Arial','fontweight','bold');

% network length
figure(4)
clf
b=bar(xxA,yyA','stacked');
hold off;ylabel('Number of networks');xlabel('Network length (nm)')
set(gca,'xtick',5:5:30,'xticklabel',0.5:0.5:3);xlim([4 31]);box off
ax = axis;hold on
plot(ax(2)*[1,1],ax(3:4),'k','linewidth',1.5)
plot(ax(1:2),ax(4)*[1,1],'k','linewidth',1.5)
hold off;legend('Quadrilateral','Pentagonal','Hexagonal');legend boxoff
ax = gca;ax.LineWidth=1.5;set(gca,'FontSize',16,'FontName', 'Arial','fontweight','bold');