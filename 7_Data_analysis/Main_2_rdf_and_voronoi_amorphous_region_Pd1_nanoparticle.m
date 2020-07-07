%% Main_2_rdf_and_voronoi_amorphous_region
% calculate the radial distribution functions and Voronoi analysis
% for the amorphous region in the sample
clear
clc
addpath('./src')

% read in files: finalized atomic coordinates in Angstrom and normalized
% BOO
data=importdata('./Input/Pd1_nanoparticle/BOO_normalized.mat');
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

% set the parameters: pixel size, step size and the range for rdf
res=0.3216;
step=0.1/res;
cutoff=20/res;

% number of cores to use for parallel computing
parpool_size=12;

model=modelamorphous/res; % use amorphus atoms, convert to unit of pixels

% calculate the alpha shape of the nanoparticle
submodel=model-repmat(min(model,[],2),[1,size(model,2)])+ones(size(model));
submodel=double(submodel);
shp = alphaShape(submodel(1,:)',submodel(2,:)',submodel(3,:)',ceil(4/res));

% initialize a spherical index for later intersection calculation
[vMat,~]=spheretri(2000);
nsample=size(vMat,1);

%% perform the main rdf calculation
xx=0:step:cutoff;
gg=xx*0;
nn=gg;
if parpool_size~=0
pjob = gcp('nocreate');
if isempty(pjob)
    parpool(parpool_size)
elseif pjob.NumWorkers ~= parpool_size
    delete(pjob)
    parpool(parpool_size)
end
end

parfor i=1:size(submodel,2)
    i
    g=xx*0;
    n=g;
    for j=1:size(submodel,2)
        dis=submodel(:,i)-submodel(:,j);
        dis=norm(dis,2);
        if dis<cutoff
            ind=ceil(dis/step+0.01);
            g(ind)=g(ind)+1;
        end
    end
    for j=xx
        spoints=vMat*(j+step/2)+repmat(submodel(:,i)',[size(vMat,1),1]);
        in = inShape(shp,spoints(:,1),spoints(:,2),spoints(:,3));
        ind=round(j/step)+1;
        n(ind)=sum(in)/nsample*(j+step/2)^2;
    end
    gg=gg+g;
    nn=nn+n;
end
rdfnorm=gg./nn;
rdfnorm=rdfnorm/rdfnorm(end-2);
g_r=imgaussfilt(rdfnorm(2:end-1),1.5);
r=res*xx(2:end-1)+0.05;
% plot(r,g_r,'-','Color',[0 0.5 0],'linewidth',2)
save('./Output/Pd1_nanoparticle/rdf_amorphous_region.mat','r','g_r')
%% perform Voronoi analysis on amorphous region

% move amorphous atoms to the beginning of the array
model=double([modelamorphous,modelcrystal])';
% indices for amorphous atoms
amorInd=1:size(modelamorphous,2);

% Voronoi regulations
areacutoff=0.01; % faects area should be larger than 1% of polygon total surface area
bondlengthcutoff=3.9; % Voronoi neighbor should be within the 1st nerest shell distance (1st valley in RDF)

% perform the voronoi calcualtion
defaultFaceColor  = [0.6875 0.8750 0.8984];
[V,R] = voronoin(model);

DT = delaunayTriangulation(model(:,1),model(:,2),model(:,3));
ConnectivityList=DT.ConnectivityList;

vor=[]; vorid=[];   vornk=[];   vorarea=[]; neigh=[];
for tid=1:size(model,1)
tid
XR10 = V(R{tid},:);
if sum(isinf(XR10),[1 2 3])
    ind=sum(ConnectivityList==tid,2)~=0;
    ind=ConnectivityList(ind,:);
    ind=unique(ind(:));
    dis=vecnorm(model(ind,:)-model(tid,:),2,2);
    ind=ind(dis<bondlengthcutoff);
    neigh{tid}=model(ind,:);
    continue
end
K = convhull(XR10);

vec1=XR10(K(:,1),:)-XR10(K(:,2),:);
vec2=XR10(K(:,1),:)-XR10(K(:,3),:);
nk=cross(vec1,vec2);
nk=nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);

ind=1:size(nk,1);
num=0;
faceid=[];
facevor=[];
while ~isempty(ind)
    flag=sum(nk(ind,:).*nk(ind(1),:),2);
    faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10(K(ind,1),:)-XR10(K(ind(1),1),:)).*nk(ind(1),:),2))<1e-5;
    tempid=K(ind(faceidtemp),:);
    num=num+1;
    faceid{num}=unique(tempid(:));
    
    % sort vertices of each face to follow clockwise or counterclockwise order
    % for plotting
    center=mean(XR10(faceid{num},:),1);
    pol=XR10(faceid{num},:)-center;
    pol=pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
    npol=size(pol,1);
    Y=dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(ind(1),:),[npol,1])');
    
    D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
    D(D<0)=360+D(D<0);
    [~,sid]=sort(D);
    faceid{num}=faceid{num}(sid);
    ind(faceidtemp)=[];
end
facenk=[];
facearea=[];
for i=1:size(faceid,2)
    % calculate surface normal
    vec1=XR10(faceid{i}(1),:)-XR10(faceid{i}(2),:);
    vec2=XR10(faceid{i}(1),:)-XR10(faceid{i}(3),:);
    nk=cross(vec1,vec2);
    facenk(i,:)=nk/sqrt(sum(nk.^2));
    
    % calculate face area
    vec=XR10(faceid{i}(2:end),:)-XR10(faceid{i}(1),:);
    facearea(i)=0.5*sum(vecnorm(cross(vec(1:end-1,:),vec(2:end,:)),2,2));
    facevor{i}=XR10(faceid{i},:);
end
vorid{tid}=faceid;
vor{tid}=facevor;
vornk{tid}=facenk;
vorarea{tid}=facearea;

% remove face with small area
faceremove = find(facearea < areacutoff*sum(facearea));
ind=sum(ConnectivityList==tid,2)~=0;
ind=ConnectivityList(ind,:);
ind=unique(ind(:));
for i=faceremove
    vec=XR10(faceid{i}(1),:)-model(tid,:);
    vec=sum(vec.*facenk(i,:))*facenk(i,:)*2;
    pos=[model(tid,:)+vec;model(tid,:)-vec];
    dis1=vecnorm(model(ind,:)-pos(1,:),2,2);
    dis2=vecnorm(model(ind,:)-pos(2,:),2,2);
    atomremove=find(dis1<1e-5 | dis2<1e-5);
    ind(atomremove)=[];
end
dis=vecnorm(model(ind,:)-model(tid,:),2,2);
ind=ind(dis<bondlengthcutoff); % comment this line if don't want bond length regulation to the voronoi
neigh{tid}=model(ind,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% re-calculate voronoi after regulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(ind)<5
    vor{tid}=[];
    continue
end
tidsub=find(ind==tid);
[Vsub,Rsub] = voronoin(model(ind,:));
XR10sub = Vsub(Rsub{tidsub},:);
if sum(isinf(XR10sub),[1 2 3])
    vor{tid}=[];
    continue
end
Ksub = convhull(XR10sub);

vec1=XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,2),:);
vec2=XR10sub(Ksub(:,1),:)-XR10sub(Ksub(:,3),:);
nk=cross(vec1,vec2);
nk=nk./repmat(sqrt(sum(nk.^2,2)),[1 3]);

indsub=1:size(nk,1);
num=0;
faceid=[];
facevor=[];
while ~isempty(indsub)
    flag=sum(nk(indsub,:).*nk(indsub(1),:),2);
    faceidtemp = abs(flag)>1-1e-5 & abs(sum((XR10sub(Ksub(indsub,1),:)-XR10sub(Ksub(indsub(1),1),:)).*nk(indsub(1),:),2))<1e-5;
    tempid=Ksub(indsub(faceidtemp),:);
    num=num+1;
    faceid{num}=unique(tempid(:));
    
    % sort vertices of each face to follow clockwise or counterclockwise order
    % for plotting
    center=mean(XR10sub(faceid{num},:),1);
    pol=XR10sub(faceid{num},:)-center;
    pol=pol./repmat(sqrt(sum(pol.^2,2)),[1 3]);
    npol=size(pol,1);
    Y=dot(cross(repmat(pol(1,:),[npol,1]),pol)',repmat(nk(indsub(1),:),[npol,1])');
    
    D = atan2d(Y,dot(repmat(pol(1,:),[npol,1])',pol'));
    D(D<0)=360+D(D<0);
    [~,sid]=sort(D);
    faceid{num}=faceid{num}(sid);
    indsub(faceidtemp)=[];
    facevor{num}=XR10sub(faceid{num},:);
end
vor{tid}=facevor;
neigh{tid}=model(ind,:);
end
% calculate indices
edgelengthcutoff=0.00; % edge length must be larger than this value

indlist=zeros(1,6); % voronoi index list
Nindlist=0; % counts of voronoi index
VoronoiID=[];
VoronoiID{1}=0; % voronoi cell IDs in each index
badID=[];
facecounts=0; % every voronoi cell has the same weight, not every face
for tid=amorInd % only count amorphous region
    facevor=vor{tid};
    if isempty(facevor)
        continue
    end
    edgelength=[];
    nedge=0;
    totedgelen=0;
    for i=1:size(facevor,2)
        ntemp=size(facevor{i},1);
        nedge=nedge+ntemp;
        edgelength{i}=vecnorm(facevor{i}-facevor{i}([ntemp 1:ntemp-1],:),2,2);
        if sum(edgelength{i}>6)
            badID=[badID tid];
            continue
        end
        totedgelen=totedgelen+sum(edgelength{i});
    end
    ind=zeros(1,6);
    for i=1:size(facevor,2)
        n=sum(edgelength{i}>=totedgelen/nedge*edgelengthcutoff);
        if n<=6
            ind(n)=ind(n)+1;
        end
    end
    facecounts=facecounts+ind/size(neigh{tid},1);
    
    temp=indlist==ind;
    id=find(sum(temp,2)==6);
    if ~isempty(id)
        Nindlist(id)=Nindlist(id)+1;
        VoronoiID{id}=[VoronoiID{id} tid];
    else
        indlist=[indlist;ind];
        Nindlist=[Nindlist,1];
        VoronoiID{end+1}=[tid];
    end
end

% plot and save results
NNN=20;
[Nindlistsort,ind]=sort(Nindlist,'descend');
tot=sum(Nindlistsort);
top10=indlist(ind(1:NNN),3:end);
Ntop10=Nindlistsort(1:NNN)/tot;

topOrder=1:NNN;
[~,ind]=sort(topOrder,'ascend');

% plot most populated Voronoi cells
figure(1)
bar(Ntop10(ind))
ylabel('Fraction')
temp=['<'*ones(NNN,1) num2str(top10(ind,:)) '>'*ones(NNN,1)];
label=cell(NNN,1);
for i=1:NNN
    label{i}=temp(i,:);
end

label = cellfun(@(x) strrep(x,'  ',','), label,'UniformOutput',false);
set(gca,'xtick',1:20,'xticklabel',label)
xtickangle(45)
disp('Top 10 most populated indices')

% 3-6 fold facets fractions
temp=indlist.*Nindlist';
figure(2) % every face has same weight
bar(sum(temp(:,3:6),1)/sum(temp(:)))
ylabel('Fraction')
save('./Output/Pd1_nanoparticle/Voronoi_amorphous_region.mat','Nindlist','indlist')