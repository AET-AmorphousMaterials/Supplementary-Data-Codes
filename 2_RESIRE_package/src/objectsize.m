function [objx,objy,objz]=objectsize(dimx,dimy,dimz,vector1,vector2,vector3,angles)
objz=dimz;
npj=size(angles,1);
dxarr=zeros(1,npj);
dyarr=zeros(1,npj);

cenx=round((dimx+1)/2);
ceny=round((dimy+1)/2);
cenz=round((dimz+1)/2);

corners=[1-cenx,1-cenx,1-cenx,1-cenx,dimx-cenx,dimx-cenx,dimx-cenx,dimx-cenx;
    1-ceny,1-ceny,dimy-ceny,dimy-ceny,1-ceny,1-ceny,dimy-ceny,dimy-ceny;
    1-cenz,dimz-cenz,1-cenz,dimz-cenz,1-cenz,dimz-cenz,1-cenz,dimz-cenz;];

for i=1:npj
    R = (MatrixQuaternionRot(vector1, angles(i,1)) * MatrixQuaternionRot(vector2, angles(i,2)) * MatrixQuaternionRot(vector3, angles(i,3)))';
    rotCorners=R*corners;
    normvec=R*[0;0;1];
    kx=normvec(1)/normvec(3);
    ky=normvec(2)/normvec(3);
    x1=((1-cenz)-rotCorners(3,:))*kx+rotCorners(1,:);
    x2=((dimz-cenz)-rotCorners(3,:))*kx+rotCorners(1,:);
    
    y1=((1-cenz)-rotCorners(3,:))*ky+rotCorners(2,:);
    y2=((dimz-cenz)-rotCorners(3,:))*ky+rotCorners(2,:);
    
    dxarr(i)=max(abs([x1,x2]));
    dyarr(i)=max(abs([y1,y2]));
end
objx=ceil(max(dxarr)*2+1);
objy=ceil(max(dyarr)*2+1);
end