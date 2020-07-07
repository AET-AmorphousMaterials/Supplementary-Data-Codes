function [obj,params,atoms] = gradient_3D_HB_all_Kmeans(para, xdata, ydata)
fprintf('\nHB gradient algorithm\n');

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;
step_sz    = xdata.step_sz;
iterClass_step=xdata.iterClass_step;

model      = single(xdata.model);
atoms       = xdata.atoms;
num_atom   = numel(atoms);
num_types  = numel(unique(atoms));

[N1,N2,N3] = size(ydata);
N_s = 2*halfWidth+1;

% atom_type1 = xdata.atom_type1;
% atom_type2 = xdata.atom_type2;
% atom_type3 = xdata.atom_type3;
% X_rot = xdata.X_rot;
% Y_rot = xdata.Y_rot;
% Z_rot = xdata.Z_rot;

fixedfa = reshape( make_fixedfa_man3d([N1,N2,N3], Res, Z_arr), [N1,N2,N3] );
model = model/Res;
%
ydata = max(real( my_ifft( my_fft(ydata) ./ fixedfa )),0);


% l2_type1 =  xdata.l2_type1;
% l2_type2 =  xdata.l2_type2;
% l2_type3 =  xdata.l2_type3;

% h1 = para(1,1);
% h2 = para(1,2);
% h3 = para(1,3);
%
% b1 = para(2,1);
% b2 = para(2,2);
% b3 = para(2,3);

X_cen = single( reshape(model(1,:),[1,1,1,num_atom]) );
Y_cen = single( reshape(model(2,:),[1,1,1,num_atom]) );
Z_cen = single( reshape(model(3,:),[1,1,1,num_atom]) );
X_round = round(X_cen);
Y_round = round(Y_cen);
Z_round = round(Z_cen);

[X_crop,Y_crop,Z_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth, -halfWidth:halfWidth);
X_crop = (single(X_crop));
Y_crop = (single(Y_crop));
Z_crop = (single(Z_crop));

[xx,yy,zz] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth, -halfWidth:halfWidth);
R4 = xx.^4 + yy.^4 + zz.^4;
sum_R4 = sum(R4(:));
clear xx yy zz R4;

h = zeros(1,1,1,num_atom, 'single');
b = zeros(1,1,1,num_atom, 'single');
if size(para,2)==num_types
    for k=1:num_types
        h(atoms==k) = para(1,k);
        b(atoms==k) = para(2,k);
    end
elseif size(para,2)==num_atom
    h(:) = para(1,:);
    b(:) = para(2,:);
else
    display('error')
    return
end
b = (Res*pi)^2./b;
b_k = zeros(num_types,1);
for k=1:num_types
    b_k(k) = mean(b(atoms==k));
end

%num_atom_type = [nnz(atom==1),nnz(atom==2),nnz(atom==3)];
%t = step_sz/(N1^2*num_pj);
sum_pj = norm( ydata(:) );

%index = zeros(2,num_pj, num_atom);
index_set    = zeros(N_s^3,num_atom, 'single');


for k=1:num_atom
    indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
    indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);
    indz = Z_round(k) + (-halfWidth:halfWidth) + round((N3+1)/2);
    [xx,yy,zz] = ndgrid(indx, indy, indz);
    index_set( : , k) = sub2ind( [N1,N2,N3], xx(:),yy(:),zz(:) );
    %temp = false(N1,N1);temp(indx,indy)=1;
    %logical_set((i-1)*N1^2 + 1 : i*N1^2, k) = logical(temp(:));
    
end

% grad_h_set = zeros(N_s^3, num_atom, 'single');
% grad_b_set = zeros(N_s^3, num_atom, 'single');

for iter=1:iterations
    %tic
    obj   = zeros(N1,N2,N3,'single');
    
    X_round = round(X_cen);
    Y_round = round(Y_cen);
    Z_round = round(Z_cen);
    
    l2 =bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
        bsxfun(@plus, Y_crop, Y_round - Y_cen).^2 + ...
        bsxfun(@plus, Z_crop, Z_round - Z_cen).^2 ;
    
    l2_b        = bsxfun(@times, l2, b);
    grad_h_set  = exp(-l2_b);
    
    pj_j        = bsxfun(@times, grad_h_set, h );
    grad_b_set  = bsxfun(@times, pj_j, l2 );
    
    for k=1:num_atom
        indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
        indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);
        indz = Z_round(k) + (-halfWidth:halfWidth) + round((N3+1)/2);
        
        obj(indx,indy,indz) = obj(indx,indy,indz) + pj_j(:,:,:,k);
        %         grad_h_set(:, k)    = reshape(grad_h_set(:,:,:,k), [N_s^3,1]);
        %         grad_b_set(:, k)    = reshape(grad_b_set(:,:,:,k), [N_s^3,1]);
    end
%     obj = max( 0, real( my_ifft( my_fft(obj) .* fixedfa ) ) );
    res = obj - ydata;
    fprintf('%d.f = %.5f\n', iter, norm( res(:) ) /sum_pj );
    
    % gradient descent
    grad_h_set = reshape(grad_h_set, [N_s^3,num_atom]);
    grad_b_set = reshape(grad_b_set, [N_s^3,num_atom]);
    grad_h = reshape(sum(res(index_set).*grad_h_set),[1,1,1,num_atom]);
    grad_b = reshape(sum(res(index_set).*grad_b_set),[1,1,1,num_atom]);
    
    % uncomment the line below if you want to update b
    %fprintf('max db = %f\n',max(abs(db(:))));
    b = b + (step_sz*.5/sum_R4)  * grad_b ./ mean(h(:)).^2;
    h = h - (step_sz/(N_s^3))* grad_h;
    
    h = max(h,0); h = min(h,200); h = max(1,h);
    b = max(b,0); %b = min(b,0.12); b = max(b,0.11);
    %b(:) = mean(b(:));
    %
    for j=1:num_types, 
        b(atoms==j) = mean(b(atoms==j));
    end
    %}
    
    %{
    if mod(iter,200)==0
        H_avg = zeros(num_types,1);
        for k=1:num_types
            H_avg(k) = mean(h(atoms==k));
        end
        [atoms,H_avg] = Kmean_HB(atoms,h(:), H_avg);
        for k=1:num_types
            h(atoms==k) = H_avg(k);
        end
    end
    %}
    %
    if mod(iter,iterClass_step)==0
        H_avg = zeros(num_types,1);
        sigma = sqrt(0.5./b(:));
        h_b = h(:).* ( 2*normcdf(4,0,sigma) - 1 );
        if nnz(h_b<=0), error('error'); end        
        for k=1:num_types            
            H_avg(k) = mean(h_b(atoms==k));            
        end
        [atoms,~] = Kmean_HB(atoms,h_b, H_avg);
%         for k=1:num_types
%             h(atoms==k) = mean(h(atoms==k));
%         end
    end
    %}
    %{
    if mod(iter,100)==0
        sigma = sqrt(0.5./b(:));
        h_b = h(:).* ( 2*normcdf(halfWidth,0,sigma) - 1 );        
        [~,sort_h] = sort(h_b);
%         atoms(sort_h(1:8859))       = 1;
%         atoms(sort_h(8860:17305))   = 2;
%         atoms(sort_h(17306:end))    = 3;
        atoms(sort_h(1:1177))       = 1;
        atoms(sort_h(1178:2571))   = 2;
        atoms(sort_h(2572:end))    = 3;
    end
    %}
    %{
    if iter>800
        for k=1:num_types
            h(atoms==k) = mean(h(atoms==k));
        end
    end
    %}
    %toc
end
%%
obj = max( 0, real( my_ifft( my_fft(obj) .* fixedfa ) ) );
params = [h(:)'; (Res*pi)^2 ./ b(:)'];

end

function [typesori,avg_H] = Kmean_HB(typesori,Hori, avg_H)


Num_species = length(avg_H);

Hcutoff=mean(Hori);

for iter1=1:1000
ind=Hori<Hcutoff;
H=Hori(ind);
types=typesori(ind);
Num_atom = length(H);
for iter=1:200
    
    obj = 0;
    old_types = types;
    
    for n = 1:Num_atom
        dist_mat = ( H(n) - avg_H ).^2;
        [dist2, idx] = min(dist_mat);
        
        types(n) = idx;
        obj = obj + dist2;
    end
    
    for k=1:Num_species
        avg_H(k) = sum(H(types==k))/sum(types==k);
    end
    
    if mod(iter,5)==1
        fprintf('%d.obj = %f\n',iter,obj/Num_atom);
    end
    
    if ~any(old_types-types), break; end
end
Hcutoff=Hcutoff-(sum(types==1)-1731)/6000*2;
if abs(sum(types==1)-1731)<=0;break;end
end
typesori(ind)=types;
typesori(~ind)=2;
end


