function [projection,interR] = calculate3Dprojection_rec_fast(modelK,phi,theta,psi,vec1,vec2,vec3,interR)

%get dimensions and centers
[dimx, dimy, dimz] = size(modelK);

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
ncz = round((dimz+1)/2);

[Y, X, Z] = meshgrid(((1:dimy) - ncy), (1:dimx) - ncx, 0);

%calculate rotation matrix
if nargin>4
    R = (MatrixQuaternionRot(vec1, phi) * MatrixQuaternionRot(vec2, theta) * MatrixQuaternionRot(vec3, psi))';
else
    R = [ cosd(psi)*cosd(theta)*cosd(phi)-sind(psi)*sind(phi) ,cosd(psi)*cosd(theta)*sind(phi)+sind(psi)*cosd(phi)   ,    -cosd(psi)*sind(theta);
         -sind(psi)*cosd(theta)*cosd(phi)-cosd(psi)*sind(phi), -sind(psi)*cosd(theta)*sind(phi)+cosd(psi)*cosd(phi) ,   sind(psi)*sind(theta)  ;
          sind(theta)*cosd(phi)                               , sind(theta)*sind(phi)                                ,              cosd(theta)];

end

R(3,:)=[];
scale = (dimy)/(dimx);
scale_mat = [1, 1/scale; scale,1; 1, 1/scale]; 

%rotate coordinates
rotkCoords = (R'.*scale_mat)*[X(:)';Y(:)'];
rotKx = rotkCoords(1,:);
rotKy = rotkCoords(2,:);
rotKz = rotkCoords(3,:)*dimz/dimx;

%reshape for interpolation
rotKx = reshape(rotKx,size(X));
rotKy = reshape(rotKy,size(Y));
rotKz = reshape(rotKz,size(Z));

if isempty(interR)
    P=[2 1 3];
    modelK = permute(modelK, P);
    interR=griddedInterpolant(modelK,'linear','none');
end

%calculate points on central slice
pjK=interR(rotKy+ncy,rotKx+ncx,rotKz+ncz);
%remove any nan from interpolation
pjK(isnan(pjK))=0;

%take IFFT to obtain projection
projection = real(my_ifft(pjK(:,:)));
end