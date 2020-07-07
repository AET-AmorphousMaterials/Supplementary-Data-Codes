function obj = interp_pj_realspace(obj)
tic
projections = obj.InputProjections;
angles = obj.InputAngles;
Num_pj = obj.NumProjs;
dimx = obj.Dim1;
dimy = obj.Dim2;
obj_dimx = obj.obj_dimx;
obj_dimy = obj.obj_dimy;
obj_dimz = obj.obj_dimz;

vec1 = obj.vector1;
vec2 = obj.vector2;
vec3 = obj.vector3;

dtype = obj.dtype;
sigma = obj.sigma;
n1_oversampled = obj.n1_oversampled;
n2_oversampled = obj.n2_oversampled;
n3_oversampled = obj.n3_oversampled;

if sigma~=0
    [Y_big,X_big,Z_big] = meshgrid(1:n2_oversampled,1:n1_oversampled,1:n3_oversampled,dtype);
    x_cen = floor(obj.n2_oversampled,dtype/2);
    y_cen = floor(obj.n1_oversampled,dtype/2);
    z_cen = floor(obj.n3_oversampled,dtype/2);
    kernel = ((X_big-x_cen).^2 + (Y_big-y_cen).^2 + (Z_big-z_cen).^2)/z_cen^2;
    obj.kernel = exp(-kernel*sigma^2);
end
clear X_big Y_big Z_big
ncx = round((dimx+1)/2); 
ncy = round((dimy+1)/2); 
k1 = cast((-1*ceil((obj_dimx-1)/2):1:floor((obj_dimx-1)/2)),dtype );
k2 = cast((-1*ceil((obj_dimy-1)/2):1:floor((obj_dimy-1)/2)),dtype );
k3 = cast((-1*ceil((obj_dimz-1)/2):1:floor((obj_dimz-1)/2)),dtype );
[XX,YY,ZZ] = ndgrid(k1,k2,k3);
XX = XX(:)';
YY = YY(:)';
ZZ = ZZ(:)';
rot_pjs = zeros(obj_dimx,obj_dimy,obj_dimz,Num_pj,dtype);

for k = 1:Num_pj
    phi   = angles(k,1);
    theta = angles(k,2);
    psi   = angles(k,3);
    pj    = double(projections(:,:,k));    
    
    R1 = MatrixQuaternionRot(vec1,phi);
    R2 = MatrixQuaternionRot(vec2,theta);
    R3 = MatrixQuaternionRot(vec3,psi);
    R =(R1*R2*R3)';
    
    rotCoords = R(1:2,:)*[XX; YY; ZZ];
    rot_x  = double(rotCoords(1,:));
    rot_y  = double(rotCoords(2,:));
    rot_pj = splinterp2(pj, rot_y+ncy, rot_x+ncx);   
    rot_pj = reshape(rot_pj, [obj_dimx, obj_dimy, obj_dimz]);
    rot_pjs(:,:,:,k) = rot_pj;        
end
obj.sum_rot_pjs = sum(rot_pjs,4);

timeTakenToInterp = toc;
timeTakenToInterp = round(10*timeTakenToInterp)./10;
fprintf('RESIRE: projections interpolated in %.12g seconds.\n\n',timeTakenToInterp);

end