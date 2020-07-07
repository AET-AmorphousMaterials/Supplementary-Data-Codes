function Projs = calcLinearPjs(para,xdata)
%fprintf('calcLinearPjs called\n');

Z_arr	= xdata.Z_arr;
Res     = xdata.Res;
N1 = xdata.volSize;
halfWidth = xdata.halfWidth;

atom_type1 = xdata.atom_type1;
atom_type2 = xdata.atom_type2;
atom_type3 = xdata.atom_type3;
X_rot = xdata.X_rot;
Y_rot = xdata.Y_rot;
Z_rot = xdata.Z_rot;
num_pj = xdata.num_pj;
num_atom = xdata.num_atom;
atom = xdata.atoms;

% l2_type1 =  xdata.l2_type1;
% l2_type2 =  xdata.l2_type2;
% l2_type3 =  xdata.l2_type3;

htAr = para(1,:);
bfAr = para(2,:);

h1 = para(1,1);
h2 = para(1,2);
h3 = para(1,3);

b1 = para(2,1);
b2 = para(2,2);
b3 = para(2,3);


[Y_crop,X_crop] = meshgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);      
Z_crop = -halfWidth:halfWidth;

Y_crop = single(Y_crop(:));
X_crop = single(X_crop(:));
Z_crop = single(Z_crop(:));

Projs = zeros(N1,N1,num_pj);
%tic
for i=1:num_pj
    for j=1:3
        atom_type_j = atom==j;
        h_j         = para(1,j);
        b_j         = para(2,j);
        %c_j = para(3,j);
        
        X_cen = X_rot(i,atom_type_j);
        Y_cen = Y_rot(i,atom_type_j);
        Z_cen = Z_rot(i,atom_type_j);
        X_round = round(X_cen);
        Y_round = round(Y_cen);
        Z_round = round(Z_cen);
        
        l2_xy = bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
                bsxfun(@plus, Y_crop, Y_round - Y_cen).^2;
        
        coef_z  = sum( exp( -bsxfun(@plus, Z_crop, Z_round - Z_cen).^2 /b_j )  );
        
        %size(l2_xy)
        %size(coef_z)
        
        pj_j = h_j* bsxfun(@times, exp(-l2_xy/b_j ), coef_z ) ;
        pj_j = reshape(pj_j, [2*halfWidth+1, 2*halfWidth+1, length(X_cen)]);
        
        for k=1:length(X_cen)
            indx = X_round(k) + (-halfWidth:halfWidth) + N1/2 + 1;
            indy = Y_round(k) + (-halfWidth:halfWidth) + N1/2 + 1;
            %         size(indx)
            %         size(pj_i)
            %size(Projs(indx,indy,k))
            %size(pj_i(:,:,k))
            %[min(indx),max(indx)]
            %[min(indy),max(indy)]
            Projs(indx,indy,i) = Projs(indx,indy,i) + pj_j(:,:,k);
        end
    end
end
%toc

% fixedfa_big = make_fixedfa_man(N1, Res, Z_arr);
% 
% Projs = cal_Vol_Projs(...
%     model,atoms,[htAr 0],[bfAr 1],N1, Res, CropHalfWidth, angles',fixedfa_big);




end