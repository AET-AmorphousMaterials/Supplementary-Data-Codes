%%  alignByNormXCorr %%

%% Align projections using normalized cross correlation
%%inputs:
%%  ref_img           - reference image; corresponds to the input projection
%%  img               - image to compare; corresponds to the calculated backprojection
%%  smooth_factor     - smoothing parameter for smooth3D. The backprojection is smoothed and thresholded to form a template
%%  threshhold_factor - value between 0 and 1 representing the fraction of maximum intensity in backprojection to use as a threshhold

%%outputs:
%%  XC - the maximum value of cross correlation
%%  new_x_center - the x (dimension 1) location of the best template match in ref_img
%%  new_y_center - the y (dimension 2) location of the best template match in ref_img

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function [metrix, new_x_center, new_y_center] = alignByNormXCorrSubpixel(ref_img,img,res)
if res>=1
    num=1;
else
    num=ceil(-log(res)/log(2));
end
t_sX=0;
t_sY=0;
metrix=0;
template=img;
for i=1:num
    for sX=t_sX+[-2^-i,0,2^-i]
        for sY=t_sY+[-2^-i,0,2^-i]
            sref_img=real(My_FourierShift(ref_img,sX,sY));
            xc = normxcorr2(template,sref_img);
            XC = max(xc(:));
            if XC>metrix
                metrix=XC;
                [xpeak, ypeak]   = find(xc==max(xc(:)));
                x_offSet = xpeak-size(template,1);
                y_offSet = ypeak-size(template,2);
                new_x_center = x_offSet + round((size(template,1)+1)/2)-sX;
                new_y_center = y_offSet + round((size(template,2)+1)/2)-sY;
                tsX=sX;
                tsY=sY;
            end
        end
    end
    t_sX=t_sX+tsX;
    t_sY=t_sY+tsY;
end
end