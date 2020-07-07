function Fcalc = BHfactorFit(para,xdata)
fprintf('call BHfactoFit\n');

model=xdata.model;
atoms=xdata.atoms;
fa_avg=xdata.fa_avg;
kx=xdata.kx;
ky=xdata.ky;
kz=xdata.kz;

bf = model(1,:);
ht = model(1,:); 
for i = 1:numel(unique(atoms))
    bf(atoms==i) = para(1,i); 
    ht(atoms==i) = para(2,i); 
end

q2 = kx.^2 + ky.^2 + kz.^2;
ra.rx = model(1,:); ra.ry = model(2,:); ra.rz = model(3,:);

L = numel(kx);

nk = length(kx); 
b_vals =  para(1,:);
index1 = atoms==1; len_r1 = sum(index1(:));
index2 = atoms==2; len_r2 = sum(index2(:));
index3 = atoms==3; len_r3 = sum(index3(:));

xj1 = ra.rx(index1); yj1 = ra.ry(index1); zj1 = ra.rz(index1);
xj2 = ra.rx(index2); yj2 = ra.ry(index2); zj2 = ra.rz(index2);
xj3 = ra.rx(index3); yj3 = ra.ry(index3); zj3 = ra.rz(index3);

sk = 2*pi*kx; %sk = sk(:);
tk = 2*pi*ky; %tk = tk(:);
uk = 2*pi*kz; %uk = uk(:);

h1 = ht(index1); %f1 = f1(:);
h2 = ht(index2); %f2 = f2(:);
h3 = ht(index3); %f2 = f2(:);

g1 = fa_avg.* exp(-b_vals(1)*q2); g1 = g1(:);
g2 = fa_avg.* exp(-b_vals(2)*q2); g2 = g2(:);
g3 = fa_avg.* exp(-b_vals(3)*q2); g3 = g3(:);

Fz1 = nufft3d3(len_r1,xj1,yj1,zj1,h1, -1,1e-6, nk,sk,tk,uk );
Fz2 = nufft3d3(len_r2,xj2,yj2,zj2,h2, -1,1e-6, nk,sk,tk,uk );
Fz3 = nufft3d3(len_r3,xj3,yj3,zj3,h3, -1,1e-6, nk,sk,tk,uk );

Fcalc = (Fz1.*g1 + Fz2.*g2 + Fz3.*g3);

% Fcalc = zeros(L,1);
% for hh = 1:L
%     Fcalc(hh) = fa_avg(hh)*sum( ht.*exp( -2*pi*1i*(kx(hh)*ra.rx+ky(hh)*ra.ry+kz(hh)*ra.rz)-bf*q2(hh) ) );
% end
% R_factor = norm(Fcalc - Fcal,'fro')/norm(Fcalc,'fro')

% k = kfactor(Fobs,Fcalc);
% Fcalc = Fcalc*k;
% Fcalc=Fcalc(:);
end

