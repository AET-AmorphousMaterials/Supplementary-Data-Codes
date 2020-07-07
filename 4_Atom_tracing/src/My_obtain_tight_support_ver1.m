% My_obtain_tight_support_ver1
% Author: Yongsoo Yang,  UCLA Physics and Astronomy
% yongsoo.ysyang@gmail.com

function curr_Supportt = My_obtain_tight_support_ver1(RECvol)

% smooth the volume
curr_RECvol = smooth3(RECvol,'b',9);

% Otsu treshold parameter (the smaller the parameter, the more generous the
% mask, i.e., larger mask)
th_dis_r_afterav = 0.9;

% Otsu threshold
im=double(curr_RECvol/max(curr_RECvol(:))*255);
[ot,outim]=otsu_thresh_3D(im);
ot=ot*max(curr_RECvol(:))/255*th_dis_r_afterav;

curr_Support = (curr_RECvol>ot) * 1;

% make the mask quite larger
se =strel3d(18);
curr_Support = imdilate(curr_Support,se);  

% make the mask quite smaller

se =strel3d(18);
curr_Supportt = imerode(curr_Support,se);  
end