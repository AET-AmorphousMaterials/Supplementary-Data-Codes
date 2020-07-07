%% Main Atom Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

% This code trace each sub area reconstruction of Ta film, then patch them
% together, which follows the scheme of sAET.

%% Step 1: trace each sub area reconstruction
addpath('./src')
for nrecon=1:65

inputDir = '../2_RESIRE_package/Output/Ta_film_sub_area_reconstruction/';
Recon_filename = [inputDir sprintf('ReconObj_sub%d.mat',nrecon)];

% load reconstructed 3D volume
Dsetvol = importdata(Recon_filename);
dimx=Dsetvol.Dim1;
dimy=Dsetvol.Dim2;
dimz=dimx;
Dsetvol=croppedOut(Dsetvol.reconstruction,[dimx,dimy,dimz]);

% Th: intensity threshold for the local maxima pixel
% local maxima with intensity less than this value will not be traced
% because they are way too weak to become actual atoms
MaxIter = 14;   CritIter = 7;   Th = 0;

% numpeak: maxmium number of peaks to trace
numpeak=100000; Res = 0.3216/3; minDist = 2.0 / Res;
SearchRad = 3;  saveInterval = 100; 

ourputstring = ['./Output/Ta_film_sub_area_tracing/Ta_film_tracing_sub_area_' num2str(nrecon)];

BoxSize0=3; %box size used for average when sorting peaks
BoxSize1=9; %box size used to find maxima
BoxSize2=9; %box size used box for fitting off-center gauss
BoxSize3=7; %used to compute visualization matrix

% upsampling the reconstruction matrix by 3*3*3 by linear interpolation
% better to run at super conputer since the size of the interpolated
% volume will be larger than 16G
xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);

xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;

xxi = xxi(3:end);
yyi = yyi(3:end);
zzi = zzi(3:end);

[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);

Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'linear',0);
clear Xi Yi Yi
FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
FinalVol_single = single(FinalVol);

% get polynomial power array
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];                
                else
                    fitCoeff(end+1,:) = [i j k 0];                
                end
            end
        end
    end
end

% get the local maxima from the reconstruction volume
se = strel3d(3);
dilatedBW = imdilate(FinalVol,se);
maxPos = find(FinalVol==dilatedBW & FinalVol>Th);
maxVals = FinalVol(maxPos);
[~,sortInd] = sort(maxVals,'descend');
maxNum = min(numpeak,length(sortInd));
maxPos = maxPos(sortInd(1:maxNum));

fprintf(1,'numpeak = %d \n',length(maxPos));

maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];  
end
clear Dsetvol dilatedBW Xi Yi Zi

% initialize the parameters
Q = 0.5;    Alpha = 1;
cropHalfSize = SearchRad;
[X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd = find(X.^2+Y.^2+Z.^2 <=(SearchRad+0.5)^2);
XYZdata.X = X(SphereInd);
XYZdata.Y = Y(SphereInd);
XYZdata.Z = Z(SphereInd);

Orders = fitCoeff(:,1:3);
PosArr = zeros(size(maxXYZ));
TotPosArr = zeros(size(maxXYZ));

exitFlagArr = zeros(1, size(maxXYZ,1));
CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

% perform the main tracing loop
for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = FinalVol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                      exitFlagArr(i) = -3;
                    else
                      TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                    end
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
    
    if mod(i,saveInterval) == 0
        parsave(sprintf('%s_result.mat',ourputstring),'PosArr',PosArr,'TotPosArr',TotPosArr,'CoeffArr',CoeffArr,'Orders',Orders,'exitFlagArr',exitFlagArr);
    end
end

parsave(sprintf('%s_result.mat',ourputstring),'PosArr',PosArr,'TotPosArr',TotPosArr,'CoeffArr',CoeffArr,'Orders',Orders,'exitFlagArr',exitFlagArr);
end

%% Step 2: patch sub area tracing to obtain final tracing of the entire volume
clear
clc
% load location information for sub area reconstructions
para=importdata('../1_Measured_data/Ta_film/Sub_area_projections_for_sAET/alignPara.mat');
cropcen=para.CropCen;
crophalf=para.CropHalf;

dims=[401 351 341];
support=zeros([401 351 341]);

dy=4; % trimed pixels around each reconstruction
BGz=20; % thickness to exclude near the top and bottom surface

unique_y=unique(cropcen(:,2));
n_unique_y=length(unique_y);
Allatom_pos=[];
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
        
        support(xx(1):xx(2),yy(1):yy(2),zz(1)+BGz:zz(2)-BGz)=1;
        
        TracingResult = importdata(['./Output/Ta_film_sub_area_tracing/Ta_film_tracing_sub_area_' num2str(ind) '_result.mat']);
        atom_pos = TracingResult.TotPosArr(TracingResult.exitFlagArr==0,:)';
        atom_pos= (atom_pos / 3) - 2;
        atm_ind= atom_pos(1,:)>=xx1(1) & atom_pos(1,:)<=xx1(2) & atom_pos(2,:)>=yy1(1) & atom_pos(2,:)<=yy1(2) ...
             & atom_pos(3,:)<=dimz-BGz*0 & atom_pos(3,:)>=BGz*0;
        
        Allatom_pos=[Allatom_pos,atom_pos(:,atm_ind)+[xx(1)-xx1(1);yy(1)-yy1(1);zz(1)-1]];
        
    end
end
% save('./Output/Ta_film_tracing_support','support')
save('./Output/Atom_tracing_all_peaks_Ta_film','Allatom_pos')