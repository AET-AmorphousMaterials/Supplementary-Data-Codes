%% RESIRE reconstruction
function obj = reconstruct(obj)
tic;
projections = obj.InputProjections;
angles      = obj.InputAngles;
vec1        = obj.vector1;
vec2        = obj.vector2;
vec3        = obj.vector3;
Num_pj      = obj.NumProjs;
step_size   = obj.step_size;
dimx        = obj.Dim1;
dimy        = obj.Dim2;
obj_dimx    = obj.obj_dimx;
obj_dimy    = obj.obj_dimy;
obj_dimz    = obj.obj_dimz;
dtype       = obj.dtype;
iterations  = obj.numIterations;

sigma       = obj.sigma;
if sigma
    kernel  = obj.kernel;
end

sum_rot_pjs = obj.sum_rot_pjs; 

if ~isempty(obj.InitialModel)
    rec = obj.InitialModel;
    rec = croppedOut(rec,[obj_dimx,obj_dimy,obj_dimz]);
else
    rec = zeros(obj_dimx,obj_dimy,obj_dimz,dtype);
end
rec_big = zeros(obj.n1_oversampled,obj.n2_oversampled,obj.n3_oversampled,dtype);
ind_V	= My_volumn_index(size(rec_big),size(rec));
rec_big( ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
projarr = zeros(dimx,dimy,Num_pj,dtype);
gradarr = zeros(obj_dimx,obj_dimy,obj_dimz,Num_pj,dtype);
dt      = (step_size/Num_pj/dimx);

fprintf('RESIRE: Reconstructing... \n\n');

if obj.monitor_R == 1
    monitorR_loopLength = obj.monitorR_loopLength;
    errR = zeros(1,floor(iterations./monitorR_loopLength));
    Rarr_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
    Rarr2_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
end

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
k1 = cast((-1*ceil((obj_dimx-1)/2):1:floor((obj_dimx-1)/2)),dtype );
k2 = cast((-1*ceil((obj_dimy-1)/2):1:floor((obj_dimy-1)/2)),dtype );
k3 = cast((-1*ceil((obj_dimz-1)/2):1:floor((obj_dimz-1)/2)),dtype );
[XX,YY,ZZ] = ndgrid(k1,k2,k3);
XX = XX(:)';
YY = YY(:)';
ZZ = ZZ(:)';

for iter=1:iterations
    recK  = double(my_fft(rec_big));    
    
    % smoothing
    if sigma
        recK = recK.*kernel;
        rec  = real(my_ifft(recK));
        rec = croppedOut(rec,[obj_dimx,obj_dimy,obj_dimz]);        
    end
    
    % compute rotated projections via FST
    [~,interR]=calculate3Dprojection_rec_fast(recK,0,0,0,vec1,vec2,vec3,[]);
    for k = 1:Num_pj
        phi   = angles(k,1);
        theta = angles(k,2);
        psi   = angles(k,3);
        [pj_cal,~] = calculate3Dprojection_rec_fast(recK,phi,theta,psi,vec1,vec2,vec3,interR);
        pj_cal = croppedOut(pj_cal, [dimx,dimy] );

        projarr(:,:,k)=pj_cal;
        R = (MatrixQuaternionRot(vec1, phi) * MatrixQuaternionRot(vec2, theta) * MatrixQuaternionRot(vec3, psi))';
        R(3,:)=[];
        rotCoords = R*[XX; YY; ZZ];
        rot_x = double(rotCoords(1,:));
        rot_y = double(rotCoords(2,:));
        rot_pj_cal = splinterp2(pj_cal, rot_y+ncy, rot_x+ncx);
        gradarr(:,:,:,k)=reshape(rot_pj_cal,[obj_dimx, obj_dimy, obj_dimz]);             
    end
    % compute gradient & apply gradient descent
    grad = sum(gradarr,4)-sum_rot_pjs;
    rec = rec - dt*grad;
    if obj.positivity
        rec = max(0,rec);
    end
    
    % compute R factor
    if obj.monitor_R == 1 && mod(iter,monitorR_loopLength) == 0
        Rarr = zeros(Num_pj,1);
        Rarr2 = zeros(Num_pj,1);
        for i=1:Num_pj
            pj = projections(:,:,i); 
            proj_i = projarr(:,:,i);
            Rarr(i) = sum(sum( abs(proj_i - pj) ))/ sum(abs(pj(:)));
            Rarr2(i) = norm( proj_i - pj, 'fro' )/ norm( pj, 'fro' );
        end
        errR1  = mean(Rarr);
        errR2 = mean(Rarr2);
        fprintf('RESIRE: Iteration %d. Rfactor=%.4f, R2factor=%.4f \n',iter, errR1, errR2);
        errR(iter./monitorR_loopLength) = errR1;
        Rarr_record(:,iter./monitorR_loopLength) = Rarr;
        Rarr2_record(:,iter./monitorR_loopLength) = Rarr2;
    else
        fprintf('RESIRE: Iteration %d \n',iter);
    end    
    
    if obj.Plot_rec
        rec1 = sum(permute(rec,[2,3,1]),3);
        rec2 = sum(permute(rec,[1,3,2]),3);
        rec3 = sum(rec,3);
        figure(1); img(rec1,'ZY',rec2,'ZX',rec3,'XY','size',[1 3]);
        drawnow();
    end
    
    rec_big( ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end
if obj.monitor_R == 1
    obj.errR            = errR;
    obj.Rarr_record     = Rarr_record;
    obj.Rarr2_record    = Rarr2_record;
end
obj.reconstruction = rec;

reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('RESIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);

end