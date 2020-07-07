%% RESIRE reconstruction
function obj = reconstruct3D2D(obj)
tic;
pjori = obj.InputProjections;
angles      = obj.InputAngles;
mask        = obj.Mask;
Num_pj      = obj.NumProjs;
step_size   = obj.step_size;
dimx        = obj.Dim1;
dimy        = obj.Dim2;
obj_dimx    = obj.obj_dimx;
obj_dimy    = obj.obj_dimy;
obj_dimz    = obj.obj_dimz;
dtype       = obj.dtype;
iterations  = obj.numIterations;

projections = permute(pjori,[1 3 2]);

if ~isempty(obj.InitialModel)
    rec = obj.InitialModel;
    rec_vec = permute(croppedOut(rec,[obj_dimx,obj_dimy,obj_dimz]),[1 3 2]);
    rec_vec = reshape(rec_vec,[obj_dimx*obj_dimz,obj_dimy]);
else
    rec_vec = zeros(obj_dimx*obj_dimz,obj_dimy,dtype);
end
dt      = (step_size/Num_pj/obj_dimz);

fprintf('RESIRE: Reconstructing... \n\n');

if obj.monitor_R == 1
    monitorR_loopLength = obj.monitorR_loopLength;
    errR = zeros(1,floor(iterations./monitorR_loopLength));
    Rarr_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
    Rarr2_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
end

% ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
% ncz = round((dimz+1)/2);

ncx_ext = ceil((obj_dimx+1)/2); 
% ncy_ext = ceil((obj_dimy+1)/2); 
ncz_ext = ceil((obj_dimz+1)/2); 

[ZZ,XX] = meshgrid( (1:1:obj_dimz) - ncz_ext, (1:1:obj_dimx) - ncx_ext );
XX = XX(:)'; ZZ = ZZ(:)';
pt_o_ratio = 4;
num_pts = size(XX,2) ;
vec_ind = ( 1:num_pts )';
rot_mat_k = cell(Num_pj,1);
for k = 1:Num_pj
    theta = angles(k,2);
    
    % compute rotation matrix R of theta
    R = [cosd(theta), -sind(theta);
         sind(theta), cosd(theta) ];   
    rot_mat_k{k} = sparse(dimx, num_pts);
    rot_x_o = ( R(1,:) * [XX;ZZ] )' + ncx;
    for s=1:pt_o_ratio
        [s1,s2] = ind2sub( [2,2], s );
        
        rot_x_shift = R(1,:) * [ (-1)^s1 ; (-1)^s2 ] *0.25;
        rot_x = rot_x_o + rot_x_shift;

        x_foor = floor(rot_x);        
        goodInd = x_foor>=1 & x_foor<dimx;
        vec_goodInd1 = vec_ind(goodInd);
        x1 = x_foor(goodInd);
        b1 = x1 + 1 - rot_x(goodInd);
        
        goodInd = x_foor==0;
        vec_goodInd2 = vec_ind(goodInd);
        x2 = x_foor(goodInd)+1;  %x2=1
        b2 = rot_x(goodInd)   ;        
        
        goodInd = x_foor==dimx;
        vec_goodInd3 = vec_ind(goodInd);
        x3 = x_foor(goodInd);
        b3 = 1+dimx - rot_x(goodInd) ;
        
        masterSub = [ [x1;x1+1; x2;x3], ...
            [vec_goodInd1;vec_goodInd1;vec_goodInd2;vec_goodInd3] ] ;
        masterVal = [b1;1-b1; b2;b3];
        
        rot_mat_k{k} = rot_mat_k{k} + accumarray(masterSub, double(masterVal), [dimx,num_pts],[],[], true );
    end
end
A = cell2mat(rot_mat_k)/pt_o_ratio;
clear Proj_op rot_mat_k rot_x  rotCoords vec_ind vec_goodInd ...
    x1 x2 z1 z2 b1 b2 masterSub masterVal XX ZZ goodint rot_x_o

timeTakenToInterp = toc;
timeTakenToInterp = round(10*timeTakenToInterp)./10;
fprintf('RESIRE: Radon transform matrix calculated in %.12g seconds.\n\n',timeTakenToInterp);
%%
tic
% initialize object
size_rec    = size(rec_vec);
grad        = zeros( size_rec ,dtype);
pj_cals     = zeros( dimx* Num_pj, dimy);
pj_cals_tmp = zeros(size(pj_cals));
projections = reshape( projections, [dimx*Num_pj, dimy] ) ;
mask = permute(mask, [1,3,2]);
mask     = reshape( mask, [dimx*Num_pj, dimy] ) ;

Avg = gallery('tridiag', dimy, 1/8,6/8,1/8);
Avg(1,1)    =6/7; Avg(2,1)      =1/7;
Avg(end,end)=6/7; Avg(end-1,end)=1/7;


for iter=1:iterations
    % forward projections via Radon transform, 
    % i.e. compute pj_cals = Au    
    for l=1:dimy
        pj_1d = double( rec_vec(:,l) );
        pj_cals_tmp(:, l) =  A* pj_1d;
    end
    %
    pj_cals(:,1)   = (6/7)*pj_cals_tmp(:,1)   + (1/7)*pj_cals_tmp(:,2);
    pj_cals(:,end) = (6/7)*pj_cals_tmp(:,end) + (1/7)*pj_cals_tmp(:,end-1);
    for l=2:dimy-1
        pj_cals(:,l) = (1/8)* pj_cals_tmp(:,l-1) + (6/8)* pj_cals_tmp(:,l) + ...
                       (1/8)* pj_cals_tmp(:,l+1);
    end
    
    % back projections:
    % i.e compute residual = Au-b, and grad = A^T(Au-b)    
    residual = (pj_cals - projections).*mask;
    
    % compute R factor
    if obj.monitor_R == 1 && mod(iter,monitorR_loopLength) == 0
        Rarr = zeros(Num_pj,1);
        Rarr2 = zeros(Num_pj,1);
        resi = reshape(residual,[dimx,Num_pj,dimy]);
        for i=1:Num_pj
            pj = pjori(:,:,i);
            resi_i=resi(:,i,:);
            Rarr(i) = sum(abs(resi_i(:)))/ sum(abs(pj(:)));
            Rarr2(i) = norm( resi_i(:), 'fro' )/ norm( pj(:), 'fro' );
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

    residual =  residual*Avg';
    for l=1:dimy
        res_l = residual(:,l) ;
        grad(:, l) =   0.1*grad(:, l) + 0.9*(res_l'*A)'  ;
    end
    rec_vec = rec_vec - (dt)*grad;
    if obj.positivity
        rec_vec = max(0,rec_vec);
    end

    if obj.Plot_rec && mod(iter,monitorR_loopLength) == 0
        rec = reshape(rec_vec, [obj_dimx, obj_dimz, obj_dimy]) ;
        rec_crop = rec; %croppedOut(rec, [dimx, obj_dimz, dimy] );
        pjYZ = sum(permute(rec_crop,[3,2,1]),3);
        pjXY = sum(permute(rec_crop,[1,3,2]),3);
        pjXZ = sum(rec_crop,3);
        figure(1); img(pjXZ,'proj XZ',pjYZ,'proj YZ',pjXY,'proj XY','size',[1 3]);
        drawnow();
    end
end
if obj.monitor_R == 1
    obj.errR            = errR;
    obj.Rarr_record     = Rarr_record;
    obj.Rarr2_record    = Rarr2_record;
end
rec = reshape(rec_vec, [obj_dimx, obj_dimz, obj_dimy]);
obj.reconstruction = permute(rec,[1 3 2]);

reconstructionTime = toc;
reconstructionTime = round(10*reconstructionTime)./10;
fprintf('RESIRE: Reconstruction completed in %.12g seconds.\n\n',reconstructionTime);

end