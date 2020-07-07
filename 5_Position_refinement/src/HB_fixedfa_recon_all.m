function vol = HB_fixedfa_recon_all(model, atoms, h_factor,b_factor, N1, halfWidth, fixedfa)

fprintf('\n3D Gaussian model reconstruction\n');

num_atom = size(model,2);

[X_crop,Y_crop,Z_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth, -halfWidth:halfWidth);

h = zeros(1,1,1,num_atom);
b = zeros(1,1,1,num_atom);
if length(h_factor)~=num_atom
    for j=1:length(h_factor)
        h(atoms==j) = h_factor(j);
        b(atoms==j) = b_factor(j);
    end
else    
    h = reshape(h_factor, [1,1,1,num_atom] );
    b = reshape(b_factor, [1,1,1,num_atom] );
end

vol   = zeros(N1,N1,N1);

%% reconstruction
X_cen = reshape(model(1,:), [1,1,1,num_atom]);
Y_cen = reshape(model(2,:), [1,1,1,num_atom]);
Z_cen = reshape(model(3,:), [1,1,1,num_atom]);

X_round = round(X_cen);
Y_round = round(Y_cen);
Z_round = round(Z_cen);

l2 = bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
    bsxfun(@plus, Y_crop, Y_round - Y_cen).^2 + ...
    bsxfun(@plus, Z_crop, Z_round - Z_cen).^2    ;

l2_b  = bsxfun(@times, l2, b);
pj_j  = bsxfun(@times, exp(-l2_b ) , h);

for k=1:num_atom
    indx = X_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
    indy = Y_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
    indz = Z_round(k) + (-halfWidth:halfWidth) + N1/2 +1;
    
    vol(indx,indy,indz) = vol(indx,indy,indz) + pj_j(:,:,:,k);
end
    
%%
if nargin>=7
    vol = max( 0, real( my_ifft( my_fft(vol) .* fixedfa ) ) );
end

end