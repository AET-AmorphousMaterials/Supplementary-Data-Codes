%% Main_1_rdf_and_boo_entire_sample
% calculate the radial distribution functions and bond orientation order
% with all atoms inside the sample
clear
clc
addpath('./src')

path='../6_Final_coordinates/';

% read in files: finalized atomic coordinates in Angstrom
model=importdata([path, 'Final_atomic_model_Ta_film.mat']);

% set the parameters: pixel size, step size and the range for rdf
res=0.3216;
step=0.1/res;
cutoff=20/res;

% number of cores to use for parallel computing
parpool_size=12;

model=model/res; % convert to unit of pixels

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
save('./Output/Ta_film/rdf_entire_sample.mat','r','g_r')

%% perform boo calculation
% read in files: finalized atomic coordinates in Angstrom
model=importdata([path, 'Final_atomic_model_Ta_film.mat']);
res=0.3216;
cutoff=4.12; % first valley position in rdf

% calculte averaged local bond orientation order parameter
Q4bar=qnbar(4,model,cutoff);
Q6bar=qnbar(6,model,cutoff);

% reference values of crystalline lattice
Q4barb=0.0363696; Q6barb=0.510688;
Q4barf=0.190941;Q6barf= 0.574524;
Q4barh=0.09722; Q6barh=0.484762;

% normalize BOO by using fcc lattice
order=sqrt(Q4bar.^2+Q6bar.^2)/sqrt(Q4barf.^2+Q6barf.^2);
save('./Output/Ta_film/BOO_normalized.mat','model','order')

% plot BOO histogram
del=0.01;
figure(1)
hist3([Q4bar;Q6bar]','nbins',[200 100],'CDataMode','auto','FaceColor','interp','LineStyle','none');
view(0,90)
xlim([0 0.3])
ylim([0 0.65])
view(0,90)
myColorMap = jet(256);
myColorMap(1,:)=1;
colormap(myColorMap);
axis tight on
cb = colorbar();
hold on
plot3([0.209 0.209],[0 0.65],[100 100],'k-','LineWidth',2)
plot3([0 0.21],[0.65 0.65],[100 100],'k-','LineWidth',2)
x=0:0.01:0.3;
y=sqrt(0.25*(Q4barf^2+Q6barf^2)-x.^2);
plot3(x,y,90*ones(size(x)),'--r','linewidth',2)
plot(Q4barb,Q6barb,'ko','markerfacecolor','k')
text(Q4barb,Q6barb+del,{'BCC'},'VerticalAlignment','bottom','HorizontalAlignment','center')
plot(Q4barf,Q6barf,'ko','markerfacecolor','k')
text(Q4barf,Q6barf+del,100,{'FCC'},'VerticalAlignment','bottom','HorizontalAlignment','center')
plot(Q4barh,Q6barh,'ko','markerfacecolor','k')
text(Q4barh,Q6barh+del,100,{'HCP'},'VerticalAlignment','bottom','HorizontalAlignment','center')
hold off
xlabel('Q4','FontSize',12)
ylabel('Q6','FontSize',12)
set(gca,'FontSize',16,'FontName', 'Arial');
xlim([0 0.21])
ylim([0 0.65])
view(0,90)
grid off
box off
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth=2;