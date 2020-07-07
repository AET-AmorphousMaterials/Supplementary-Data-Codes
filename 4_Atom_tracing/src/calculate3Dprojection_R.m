function projection = calculate3Dprojection_R(Recon,phi,theta,psi, custom_euler_rot_vecs, imx,imy)

%get dimensions and centers
[dimx, dimy, dimz] = size(Recon);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2);

[Y, X, Z] = meshgrid((1:dimy) - ncy, (1:dimx) - ncx, (1:dimz) - ncz);

%calculate rotation matrix
R = (MatrixQuaternionRot(custom_euler_rot_vecs{1}, phi) * MatrixQuaternionRot(custom_euler_rot_vecs{2}, theta) * MatrixQuaternionRot(custom_euler_rot_vecs{3}, psi));

%[ky, kx, kz ] = meshgrid((1:dimy) - ncy, (1:dimx) - ncx, (1:dimz) - ncz);

%rotate coordinates
rotkCoords = R*[X(:)';Y(:)';Z(:)'];
rotKx = rotkCoords(1,:);
rotKy = rotkCoords(2,:);
rotKz = rotkCoords(3,:);

%reshape for interpolation
rotKx = reshape(rotKx,size(X));
rotKy = reshape(rotKy,size(Y));
rotKz = reshape(rotKz,size(Z));

%calculate points on central slice
%pjK = interp3(ky,kx,kz,double(Recon),rotKy,rotKx,rotKz,'linear');
%pjK = interp3(ky,kx,kz,Recon,rotKy,rotKx,rotKz,'cubic');
 pjK = splinterp3(double(Recon),rotKy+ncy,rotKx+ncx,rotKz+ncz);
%remove any nan from interpolation
pjK(isnan(pjK))=0;
% pjK(isinf(pjK))=0;
%take IFFT to obtain projection
projection = sum(pjK,3);

cenimx=round((imx+1)/2);
cenimy=round((imy+1)/2);
projection = projection((1:imx)-cenimx+ncx,(1:imy)-cenimy+ncy);

end