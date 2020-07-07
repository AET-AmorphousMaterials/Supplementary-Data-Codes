function [Projs,params,errR] = gradient_H_all(para, xdata, ydata)
errR=[];
fprintf('\nHB gradient algorithm\n');

Z_arr	   = xdata.Z_arr;
Res        = xdata.Res;
halfWidth  = xdata.halfWidth;
iterations = xdata.iterations;
step_sz    = xdata.step_sz;

model      = single(xdata.model);
angles     = xdata.angles;
atom       = xdata.atoms;
num_atom   = numel(atom);
num_atom_type = numel(unique(atom));

[N1,N2,num_pj] = size(ydata);
N_s = 2*halfWidth+1;

% atom_type1 = xdata.atom_type1;
% atom_type2 = xdata.atom_type2;
% atom_type3 = xdata.atom_type3;
% X_rot = xdata.X_rot;
% Y_rot = xdata.Y_rot;
% Z_rot = xdata.Z_rot;

fixedfa = reshape( make_fixedfa_man([N1,N2], Res, Z_arr), [N1,N2] );
model = model/Res;
% 
% for i=1:num_pj
%     ydata(:,:,i) = max(real( my_ifft( my_fft( ydata(:,:,i) ) ./fixedfa )),0);
% end

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

X_rot = zeros(num_pj,num_atom,'single');
Y_rot = zeros(num_pj,num_atom,'single');
Z_rot = zeros(num_pj,num_atom,'single');

for i=1:num_pj
  R1 = MatrixQuaternionRot([0 0 1],angles(i,1));  
  R2 = MatrixQuaternionRot([0 1 0],angles(i,2));
  R3 = MatrixQuaternionRot([1 0 0],angles(i,3));  
  R   = (R1*R2*R3)';
      
  rotCoords = R*model;
  X_rot(i,:) = rotCoords(1,:);
  Y_rot(i,:) = rotCoords(2,:);
  Z_rot(i,:) = rotCoords(3,:);
end


[X_crop,Y_crop] = ndgrid( -halfWidth:halfWidth, -halfWidth:halfWidth);      
Z_crop = -halfWidth:halfWidth;
X_crop = (single(X_crop));
Y_crop = (single(Y_crop));
Z_crop = (single(Z_crop));

h = zeros(1,1,num_atom, 'single');
b = zeros(1,1,num_atom, 'single');
if size(para,2)==num_atom_type
    for k=1:num_atom_type
        h(atom==k) = para(1,k);
        b(atom==k) = para(2,k);
    end
elseif size(para,2)==num_atom
    h(:) = para(1,:);
    b(:) = para(2,:);
else
    display('error')
    return
end
b = (Res*pi)^2./b;

%num_atom_type = [nnz(atom==1),nnz(atom==2),nnz(atom==3)];
% t = step_sz/(N1^2*num_pj);
% sum_pj = norm( ydata(:) );

%index = zeros(2,num_pj, num_atom);
index_set    = zeros(num_pj*N_s^2,num_atom, 'single');
%logical_set  = false(num_pj*N1^2,num_atom);
grad_h_set = zeros(N_s^2*num_pj, num_atom, 'single');
grad_b_set = zeros(N_s^2*num_pj, num_atom, 'single');

for i=1:num_pj
    X_round = round(X_rot(i,:));
    Y_round = round(Y_rot(i,:));
    
    for k=1:num_atom
%         [i,k]
        indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
        indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);
        [yy,xx] = meshgrid(indy, indx);
        index_set( (i-1)*N_s^2 + 1 : i*N_s^2, k) = sub2ind( [N1,N2,num_pj], xx(:),yy(:),i*ones(N_s^2,1) );
        %temp = false(N1,N1);temp(indx,indy)=1;
        %logical_set((i-1)*N1^2 + 1 : i*N1^2, k) = logical(temp(:));
        
    end
end

for iter=1:iterations
Projs   = zeros(N1,N2,num_pj,'single');

%tic
for i=1:num_pj
    
%   RM1 = MatrixQuaternionRot([0 0 1],angles(1,i));  
%   RM2 = MatrixQuaternionRot([0 1 0],angles(2,i));
%   RM3 = MatrixQuaternionRot([1 0 0],angles(3,i));
%   model_rot = (RM1*RM2*RM3)'*model; 
%   x_rot = model_rot(1,:);
%   y_rot = model_rot(2,:);
%   z_rot = model_rot(3,:);
  
        
        X_cen = reshape(X_rot(i,:), [1,1,num_atom]);
        Y_cen = reshape(Y_rot(i,:), [1,1,num_atom]);
        Z_cen = reshape(Z_rot(i,:), [1,1,num_atom]);
        
        X_round = round(X_cen);
        Y_round = round(Y_cen);
        Z_round = round(Z_cen);
        
        l2_xy = bsxfun(@plus, X_crop, X_round - X_cen).^2 + ...
                bsxfun(@plus, Y_crop, Y_round - Y_cen).^2;
        
        l2_z  = bsxfun(@plus, Z_crop, Z_round - Z_cen).^2    ;
        
%         size(l2_xy)
%         size(l2_z)
%         size(b)        

        l2_xy_b     = bsxfun(@times, l2_xy, b);
        l2_z_b      = bsxfun(@times, l2_z , b);        
        exp_l2_xy_b = exp(-l2_xy_b);
        exp_l2_z_b  = exp(-l2_z_b);
        
        pj_j     = bsxfun(@times, exp_l2_xy_b, sum(exp_l2_z_b) );
        %t_h = sum(sum(pj_j,1),2);
        pj_j_h   = bsxfun(@times, pj_j , h);        
        
        grad_exp = bsxfun(@times, exp_l2_xy_b, sum(l2_z.*exp_l2_z_b) );        
        bj_j     = bsxfun(@times, h, grad_exp) + pj_j_h.*l2_xy ;
        
        %pj_j = reshape(pj_j, [2*halfWidth+1, 2*halfWidth+1, length(X_cen)]);
        
        for k=1:num_atom
            indx = X_round(k) + (-halfWidth:halfWidth) + round((N1+1)/2);
            indy = Y_round(k) + (-halfWidth:halfWidth) + round((N2+1)/2);

            Projs(indx,indy,i) = Projs(indx,indy,i) + pj_j_h(:,:,k);
            %grad_h(:,:,i,k) = pj_j(:,:,k);
            grad_h_set((i-1)*N_s^2 + 1 : i*N_s^2, k) = reshape(pj_j(:,:,k), [N_s^2,1]);
            %index(:,i,k) = [X_round(k);Y_round(k)];
            grad_b_set((i-1)*N_s^2 + 1 : i*N_s^2, k) = reshape(bj_j(:,:,k), [N_s^2,1]);
        end
    
    %Projs(:,:,i) = real( my_ifft( my_fft( Projs(:,:,i)) .* fixedfa ) );    
end


for i=1:num_pj
    Projs(:,:,i) = real(( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ));
end

% k=sum(Projs(:).*ydata(:))/sum(Projs(:).^2);
% Projs=Projs*k;
% h=h*k;
% grad_b_set=grad_b_set*k;
% grad_h_set=grad_h_set*k;

res = Projs - ydata;
errR(iter)=sum( abs(Projs(:)-ydata(:)) ) / sum(abs(ydata(:)));
fprintf('%d.f = %.5f\n', iter, errR(iter) );

% for i=1:num_pj
%     res(:,:,i) = (( my_ifft( my_fft(res(:,:,i)) .* conj(fixedfa) ) ));
% end

% gradient descent
grad_b = reshape(sum(res(index_set).*grad_b_set),[1,1,num_atom]);
grad_h = reshape(sum(res(index_set).*grad_h_set),[1,1,num_atom]);

% b = b + (t/halfWidth^4) * grad_b./max(1,h.^2);
%h = h - 300*t* ( grad_h );
%b = b + (1/N_s^7/num_pj) * grad_b ./max(1,h.^2);
%[N_s^4*num_pj,max(abs(grad_h(:)))]
%h = h - (1/N_s^4/num_pj) * grad_h ;
%[N_s^2*num_pj*max(abs(t_h(:))), max(abs(grad_h(:)))]
%delta_h =  (1/N_s^4/num_pj) * grad_h ;
%delta_h =  (1/N2^2/num_pj) * grad_h;
% delta_h =  (1/num_pj/N_s^2) * grad_h ./t_h;
%[min(abs(delta_h(:))), max(abs(delta_h(:))), mean(abs(delta_h(:)))]
% h = h - delta_h;

% uncomment the line below if you want to update b
% b = b + (step_sz/(N1^2*N_s^6*num_pj))  * grad_b ; 
h = h - (step_sz/(N1^2*num_pj))        * grad_h;


% for k=1:num_atom
% %     grad_h_k=0;    
% %     for i=1:num_pj
% %         indx = index(1,i,k) + (-halfWidth:halfWidth) + N1/2 +1;
% %         indy = index(2,i,k) + (-halfWidth:halfWidth) + N1/2 +1;
% %         grad_h_k = grad_h_k + sum(sum( res(indx,indy,i) .* grad_h(:,:,i,k) ));
% %     end
% %     size(res(index_set(:,k)))
% %     size(grad_h_set(:,k))
%     grad_h_k = dot( res(index_set(:,k)), grad_h_set(:,k));
%     %grad_h_k = dot( res(logical_set(:,k)), grad_h_set(:,k));    
%     h(k) = h(k) - (t/h(k)) * grad_h_k;
%     
%     grad_b_k = dot( res(index_set(:,k)), grad_b_set(:,k) );
%     %0.00001*(t/max(h)) * grad_b_k
%     b(k) = b(k) + (t/h(k)^2/halfWidth^4) * grad_b_k;
% end

h = max(h,0);
b = max(b,0);
%fprintf('h = [%s,%s,%s]\n',h);
%toc
end
%%
% for i=1:num_pj
%     Projs(:,:,i) = max( 0, real( my_ifft( my_fft(Projs(:,:,i)) .* fixedfa ) ) );
% end
params = [h(:)'; (Res*pi)^2 ./ b(:)'];
end