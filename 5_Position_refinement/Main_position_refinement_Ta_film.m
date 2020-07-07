%% Main Position refinement 
% refinement of the atomic coordinates by minimizing the error between the
% calculated linear projections and the measured projections
clear
clc
addpath('./src/')
inputDir = './Input/Ta_film/';

% read in files: measured projections, angles and atomic position
angles=importdata([inputDir 'angles.mat']);
projections=importdata([inputDir 'projections.mat']);
model=importdata([inputDir 'atomic_model_Ta_film.mat']);

% process and select input data
pjs = [1:46];
angles = angles(pjs,:);
projections = max(projections(:,:,pjs),0);
projections=My_paddzero(projections,size(projections)+[50 50 0]);

atoms = ones(1,size(model,2));
[N1,N2,num_pj] = size(projections);
halfWidth = 4; % the cropped bondary size for each atoms
Z_arr = [73]; % atomic Z number of Ta
Res = 0.3216; % indicate the pixel size for measured projections

% initialize refinement parameters
xdata = [];
xdata.Res     = Res;
xdata.Z_arr   = Z_arr;
xdata.halfWidth = halfWidth;
xdata.atoms = atoms;
xdata.model = model;
xdata.angles = angles;

% starting H, B factor
para0 = [80;
    6];
model_refined = model;

% repeating the main refinement few times
for iiiiii=1:2
x0 = para0;
x0(1,:)=x0(1,:)/x0(1,1);

xdata.model = model_refined;
xdata.model_ori = model_refined;

% search for optimum B factor by using scanning method

% scan with large steps
num=0;
for h_scan=0
    for b_scan=-4:0.2:4
        num=num+1;
        x(:,:,num)=x0+[h_scan;b_scan];
    end
end
err_arr=zeros(1,num);
for i=1:num
    [y_pred,para0] = Cal_Bproj(x(:,:,i), xdata, projections);
    err_arr(i)=sum( abs(y_pred(:)-projections(:)) ) / sum(abs(projections(:)));
end
errRscan=[];
[errR,ind]=min(err_arr);
errRscan(end+1)=errR;

% scan with small steps
num=0;
x0=x(:,:,ind);
for h_scan=0
    for b_scan=-0.2:0.01:0.2
        num=num+1;
        x(:,:,num)=x0+[h_scan;b_scan];
    end
end
err_arr2=zeros(1,num);
for i=1:num
    [y_pred,para0] = Cal_Bproj(x(:,:,i), xdata, projections);
    err_arr2(i)=sum( abs(y_pred(:)-projections(:)) ) / sum(abs(projections(:)));
end
errRscan=[];
[errR,ind]=min(err_arr2);
errRscan(end+1)=errR;
para0=x(:,:,ind);
[y_pred,para0] = Cal_Bproj(para0, xdata, projections);
% save([inputDir 'Bscan_' num2str(iiiiii) '.mat'],'y_pred','para0','errRscan','err_arr','err_arr2')

% optimized around the above best B factor by gradient descend
xdata.step_sz    = 1;
xdata.iterations = 100;
[y_pred,para0,errR] = gradient_B(para0, xdata, projections);
% save([inputDir 'Bscan_' num2str(iiiiii) '_itr20.mat'],'y_pred','para0','errR')

% refined atomic coordinates with the best refined B factor by gradient descend
xdata.step_sz    = 1;
xdata.iterations = 1000;
[y_pred,model_arr,errR] = gradient_XYZ_man(para0, xdata, projections); model_refined = model_arr(:,:,end);
% save([inputDir 'Bxyz_' num2str(iiiiii) '_itr100.mat'],'y_pred','model_arr','errR')
end
save('./Output/Final_atomic_model_Ta_film.mat','model_refined')