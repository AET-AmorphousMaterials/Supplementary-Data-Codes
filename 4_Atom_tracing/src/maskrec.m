function maskV = maskrec(rec,lobj_dimx,hobj_dimx,lobj_dimy,hobj_dimy)
maskV=ones(size(rec));
dimx=size(maskV,1);
for i=1:size(rec,2)
    for j=1:size(rec,3)
        maskV(1:lobj_dimx,i,j)=(1:lobj_dimx)/lobj_dimx;
        maskV(hobj_dimx:end,i,j)=((dimx:-1:hobj_dimx)-hobj_dimx)/(dimx-hobj_dimx);
    end
end

nz=size(rec,3);
wz=(1:nz)-round((nz+1)/2);
wz=1-abs(wz)./nz;

for i=1:size(rec,1)
    for j=1:size(rec,2)
        maskV(i,j,:)=wz'.*squeeze(maskV(i,j,:));
    end
end
maskV(lobj_dimx:hobj_dimx,lobj_dimy:hobj_dimy,:)=1;