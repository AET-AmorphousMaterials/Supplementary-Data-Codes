%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%                     Welcome to RESIRE!                         %
%          REal Space Iterative Reconstruction Engine            %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This version incorperated extended sample recosntruction as well as 
% Radon transform for handling large samples 

addpath('./src/')
%%%%%%%%%%%%%%%%%%%%% Set RESIRE Parameters %%%%%%%%%%%%%%%%%%%%%
% See the object RESIRE_Reconstructor() for description of parameters

% please run this on super computer due to the large array size and
% oversampling ratio

% Construct all sub area regions for sAET
for i=1:65
RESIRE = RESIRE_reconstructor();
Path = '../1_Measured_data/Ta_film/';
RESIRE.filename_Projections = [Path sprintf('Sub_area_projections_for_sAET/Projections_area%d.mat',i)];
RESIRE.filename_Angles = [Path 'angles.mat'];
RESIRE.filename_Support = ''; 
RESIRE.filename_InitialModel = '';
RESIRE.filename_Results = ['./Output/Ta_film_sub_area_reconstruction/' sprintf('ReconObj_sub%d.mat',i)];

RESIRE.numIterations = 500;
RESIRE.extenedobject = true;
RESIRE.obj_dimz = 1; % object thickness, 1 - same as projection dim1; 2 - same as projection dim2; Other - specified size
RESIRE.positivity = true;
RESIRE.method = 'FST'; % 'FST' or 'Radon' method for forward projection
RESIRE.oversamplingRatio_x =1.5;
RESIRE.oversamplingRatio_y =3;
RESIRE.oversamplingRatio_z =3;
RESIRE.vector1 = [0 0 1];
RESIRE.vector2 = [0 1 0];
RESIRE.vector3 = [1 0 0];
RESIRE.step_size = 1;
RESIRE.griddingMethod = 1;
RESIRE.monitor_R = 1;
RESIRE.monitorR_loopLength = 1;
RESIRE.Plot_rec = 0;
RESIRE.dtype = 'single';
%%%%%%%%%%%%%%%%%%%%% Begin RESIRE %%%%%%%%%%%%%%%%%%%%%
RESIRE = readFiles(RESIRE);
RESIRE = CheckPrepareData(RESIRE);
RESIRE = reconstruct_control(RESIRE);
Reconstruction = RESIRE.reconstruction;
% save([Path 'Reconstruction.mat'],'Reconstruction')
SaveResults(RESIRE);
end

%% Patch sub area reconstructions to obtain final large volume
clear
clc
Path = '../1_Measured_data/Ta_film/';
% load location information for sub area reconstructions
para=importdata([Path 'Sub_area_projections_for_sAET/alignPara.mat']);
cropcen=para.CropCen;
crophalf=para.CropHalf;

dy=4; % trimed pixels around each reconstruction
Recon=zeros([401 351 341]);
unique_y=unique(cropcen(:,2));
n_unique_y=length(unique_y);
err=[];
for i=1:n_unique_y
    unique_x=cropcen(cropcen(:,2)==unique_y(i),1);
    n_unique_x=length(unique_x);
    for j=1:n_unique_x
        ind=find(cropcen(:,1)==unique_x(j) & cropcen(:,2)==unique_y(i));
        half=crophalf(ind,:);
        if i==1
            yy=[unique_y(i)-half(2),unique_y(i)+half(2)-dy];
            yy1=[1,half(2)*2+1-dy];
        elseif i==n_unique_y
            yy=[unique_y(i)-half(2)+dy,unique_y(i)+half(2)];
            yy1=[1+dy,half(2)*2+1];
        else
            yy=[unique_y(i)-half(2)+dy,unique_y(i)+half(2)-dy];
            yy1=[1+dy,half(2)*2+1-dy];
        end
        zz=[cropcen(ind,3)-crophalf(ind,3),cropcen(ind,3)+crophalf(ind,3)];
        dx=round(half(1)-unique_x(2)/2+unique_x(1)/2)-2;
        if j==1
            xx=[unique_x(j)-half(1),unique_x(j)+half(1)-dx];
            xx1=[1,half(1)*2+1-dx];
        elseif j==n_unique_x
            xx=[unique_x(j)-half(1)+dx,unique_x(j)+half(1)];
            xx1=[1+dx,half(1)*2+1];
        else
            xx=[unique_x(j)-half(1)+dx,unique_x(j)+half(1)-dx];
            xx1=[1+dx,half(1)*2+1-dx];
        end
        subrecon=importdata(['./Output/Ta_film_sub_area_reconstruction/' sprintf('ReconObj_sub%d.mat',ind)]);
%         err=[err, subrecon.errR(end)];
        [dimx,dimy,~]=size(subrecon.InputProjections);
        dimz=dimx;
        subrecon=croppedOut(subrecon.reconstruction,[dimx,dimy,dimz]);
        Recon(xx(1):xx(2),yy(1):yy(2),zz(1):zz(2))=subrecon(xx1(1):xx1(2),yy1(1):yy1(2),:);
    end
end
save('./Output/Ta_film_volume.mat','Recon')