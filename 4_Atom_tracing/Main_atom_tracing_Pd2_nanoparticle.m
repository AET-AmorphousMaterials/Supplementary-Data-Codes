%% Main Atom Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

addpath('./src')
inputDir = '';
Recon_filename = [inputDir '../3_Final_reconstruction_volume/Pd2_particle_volume.mat'];

% load reconstructed 3D volume
Dsetvol = importdata(Recon_filename);

% Th: intensity threshold for the local maxima pixel
% local maxima with intensity less than this value will not be traced
% because they are way too weak to become actual atoms
MaxIter = 14;   CritIter = 7;   Th = 1;

% numpeak: maxmium number of peaks to trace
saveInterval = 1000;    numpeak=500000; Res = 0.454/3;
minDist = 2.0 / Res;    SearchRad = 3;

ourputstring = [inputDir './Output/Atom_tracing_all_peaks_Pd2_nanoparticle'];

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

Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
clear Xi Yi Zi

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
        save(sprintf('%s.mat',ourputstring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');
    end
end
save(sprintf('%s.mat',ourputstring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');