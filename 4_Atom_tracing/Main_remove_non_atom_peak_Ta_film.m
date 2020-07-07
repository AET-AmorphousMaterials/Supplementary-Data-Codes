%% Main remove non atom peak
% Use K-mean method to remove non atom peaks from tracing results
clear
clc
addpath('./src/');
inputpath='./Input/';
saveprefix='Atom_tracing_non_atom_removed_Ta_film';

% load recontrcution volume
rec = importdata('../3_Final_reconstruction_volume/Ta_film_volume.mat');
rec=double(rec);

% load tracing results
atom_pos = importdata([inputpath 'Atom_tracing_all_peaks_Ta_film.mat']);

b1 = find(atom_pos(1,:)<1 | atom_pos(1,:)>size(rec,1)-1);
b2 = find(atom_pos(2,:)<1 | atom_pos(2,:)>size(rec,2)-1);
b3 = find(atom_pos(3,:)<1 | atom_pos(3,:)>size(rec,3)-1);
bT = union(union(b1,b2),b3);
atom_pos(:,bT) = [];

% load support to mask out the volume outside the sample boundary from the
% reconstruction
tight_support=importdata([inputpath 'tight_support_Ta_film.mat']);
se =strel3d(10);
tight_support  = imerode(tight_support ,se);  

% If there are two atoms near the boundary of sub reconstructions are close
% to each other, then keep the averaged position of the two atoms
curr_model0      = atom_pos;
count=0;
exatom=[];
i=1;
while i<=size(curr_model0,2)
    dif=curr_model0-curr_model0(:,i);
    dis=sqrt(dif(1,:).^2+dif(2,:).^2+dif(3,:).^2);
    ind=find(dis<3);
    if length(ind)>1
        count=count+length(ind)-1;
        exatom=[exatom,curr_model0(:,ind)];
        temp=mean(curr_model0(:,ind),2);
        curr_model0(:,ind)=[];
        curr_model0(:,end+1)=temp;
        i=i-1;
    end
    i=i+1;
end

% set parameters for the K-mean calculation
lnorm = 1; % error metric order
interp_type = 'cubic'; % upsampling method for the reconstruction
curr_model      = curr_model0;
labels          = ones(1,size(curr_model0,2));
Num_species     = 2; % number of different species to seperate by K-mean
Num_atom        = size(curr_model,2);
O_Ratio = 3; % up sampling ratio
halfSize = 4/O_Ratio; % number of pixels for atoms
plothalfSize = 7/O_Ratio; % number of pixels for atoms for ploting
ds = 1/O_Ratio;

% obtain global intensity histogram
[XX,YY,ZZ] = ndgrid(-halfSize: ds :halfSize, -halfSize: ds :halfSize, -halfSize: ds :halfSize);
SphereInd = find(XX.^2+YY.^2+ZZ.^2 <=(halfSize+0.5*ds)^2);

[XXp,YYp,ZZp] = ndgrid(-plothalfSize: ds :plothalfSize, -plothalfSize: ds :plothalfSize, -plothalfSize: ds :plothalfSize);
SphereInd_plot = find(XXp.^2+YYp.^2+ZZp.^2 <=(plothalfSize+0.5*ds)^2);

SPHyn = 1;
if SPHyn
    useInd      = SphereInd;
    useInd_plot = SphereInd_plot;
else
    useInd = 1:length(XX);
    useInd_plot = 1:length(XXp);
end

% generate points coordinates
YY = YY(useInd); XX = XX(useInd); ZZ = ZZ(useInd);
y_set = zeros(length(XX), Num_atom);
x_set = zeros(length(YY), Num_atom);
z_set = zeros(length(ZZ), Num_atom);

YYp = YYp(useInd_plot); XXp = XXp(useInd_plot); ZZp = ZZp(useInd_plot);
y_set_plot = zeros(length(XXp), Num_atom);
x_set_plot = zeros(length(YYp), Num_atom);
z_set_plot = zeros(length(ZZp), Num_atom);

% interpolations for points
for k=1:Num_atom
    y_set(:,k) = YY + curr_model(2,k);
    x_set(:,k) = XX + curr_model(1,k);
    z_set(:,k) = ZZ + curr_model(3,k);
    
    y_set_plot(:,k) = YYp + curr_model(2,k);
    x_set_plot(:,k) = XXp + curr_model(1,k);
    z_set_plot(:,k) = ZZp + curr_model(3,k);
end
if strcmp(interp_type,'linear')
    points = splinterp3(rec, y_set,x_set,z_set);
else
    points =    interp3(rec, y_set,x_set,z_set, interp_type);
end
points(isnan(points))=0;
points(isinf(points))=0;
points_plot = splinterp3(rec, y_set_plot, x_set_plot, z_set_plot);
points_plot(isnan(points_plot))=0;
points_plot(isinf(points_plot))=0;

% integrate intensity (for 3x3x3 and 5x5x5 voxels) for each traced peak
intensity_integ      = sum(points);
intensity_integ_max  = max(points);
intensity_integ_plot = sum(points_plot);

separate_part = 100;
L_forAver = 10;
[hist_inten,cen_integ_total]= hist(intensity_integ,separate_part);

% calculate starting reference for each species
dist_mat = zeros(Num_species,1);

[hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);

% fit histogram with two gaussian
i_guess = [max(hist_inten_plot)/10 cen_integ_total_plot(round(separate_part/6)) cen_integ_total_plot(round(separate_part/10))...
        max(hist_inten_plot) cen_integ_total_plot(round(separate_part/4*2)) cen_integ_total_plot(round(separate_part/10))];
Xdata = cen_integ_total_plot;
Ydata = hist_inten_plot;
[p, fminres, fitresult] = My_two_gaussianfit(Xdata, Ydata, i_guess);

% signle gaussian used in fitting
fun = @(p,xdata) p(1)*exp(-((xdata-p(2))/p(3)).^2);
fitresult1=fun(p(1:3),Xdata);
fitresult2=fun(p(4:6),Xdata);

ind1=find(abs(intensity_integ_plot-p(2))<p(3)/10);
ind2=find(abs(intensity_integ_plot-p(5))<p(6)/10);

avg_atom(:,1) = sum(points(:,ind1),2)/length(ind1);
avg_atom(:,2) = sum(points(:,ind2),2)/length(ind2);
avg_atom(isnan(avg_atom))=0;

y_up = round(max(hist_inten_plot)/10*12);
iter=0;

% main K-mean iteration
while true
    iter=iter+1;
    intensity_integ_1_plot = intensity_integ_plot(labels == 1);
    intensity_integ_2_plot = intensity_integ_plot(labels == 2);
    intensity_integ_3_plot = intensity_integ_plot(labels == 3);

    figure(200+2*halfSize*O_Ratio+1)
    clf
    subplot(Num_species+1,1,1);
    hist(intensity_integ_plot,separate_part);
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('OR=%d, interp=%s\nboxsize %d',O_Ratio,interp_type,2*halfSize*O_Ratio+1));
    
    subplot(Num_species+1,1,2);
    hist(intensity_integ_1_plot,cen_integ_total_plot);
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d Non atoms',sum(labels==1)));
    
    subplot(Num_species+1,1,3);
    hist(intensity_integ_2_plot,cen_integ_total_plot)
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d Ta atoms',sum(labels==2)));
      
    avg_atom1 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom2 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);

    avg_atom1(useInd) = avg_atom(:,1);
    avg_atom2(useInd) = avg_atom(:,2);
    figure(300+2*halfSize*O_Ratio+1);
    clf
    img(sum(avg_atom1,3),'atom1', sum(avg_atom2,3),'atom2');
    drawnow;
    
    old_labels = labels;
       
    obj = 0;
    for n = 1:Num_atom
        for k=1:Num_species
            dist_mat(k) = norm(points(:,n) - avg_atom(:,k), lnorm).^lnorm/sum(points(:,n).^lnorm);
            if sum(points(:,n).^lnorm)==0
                dist_mat(k)=0;
            end
        end
        [dist2, idx] = min(dist_mat);        
        labels(n) = idx;
        if isinf(dist2)
            n
        end
        obj = obj + dist2;
    end
    
    for k=1:Num_species
        avg_atom(:,k) = sum(points(:,labels==k),2)/sum(labels==k);
    end
    avg_atom(:,1)=0; % constraint non atom peak to be zero
    avg_atom(isnan(avg_atom))=0;
    fprintf('%02i. obj = %.3f\n',iter, (obj/Num_atom)^(1/lnorm) );
    
    % if there is no change in the atomic specise classification, break
    if ~any(old_labels-labels), break; end
    
end

% apply support
atomtype=zeros(1,size(curr_model,2));
for i=1:length(atomtype)  %if previous determined noatom is within the tight support, then change back to atom
        Rpos = round(curr_model(:,i));
        if tight_support(Rpos(1),Rpos(2),Rpos(3))>0
            atomtype(i) = 1;
        else
            atomtype(i) = 0;  %set atom outside support to non atoms
        end
end
atomtype=atomtype==0;
curr_model(:,atomtype)=[];
labels(atomtype)=[];

save(['./Output/' saveprefix],'curr_model','labels','intensity_integ_plot','cen_integ_total_plot','intensity_integ_max','hist_inten_plot');