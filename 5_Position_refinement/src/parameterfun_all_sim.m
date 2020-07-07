function [dd1, dd2, k] = parameterfun_all_sim(x,ra,fa,ka,L,Fobs,Fcalc,atoms)
% Yao Yang
% fit to any number species

Num_type = size(fa,1) - 1;
fa_n = fa(Num_type,:);

bf = zeros(1,numel(atoms));
ht = zeros(1,numel(atoms));
for i = 1:Num_type
bf(atoms==i) = x(1,i);
ht(atoms==i) = x(2,i);
end

ka_kx = ka.kx; ka_ky = ka.ky; ka_kz = ka.kz;
ra_rx = ra.rx; ra_ry = ra.ry; ra_rz = ra.rz;

M = size(ra_rx,2);
dbf = zeros(1,M);
dht = zeros(1,M);
s2 = ka_kx.^2+ka_ky.^2+ka_kz.^2;

parfor hh = 1:L
Fcalc(hh) = fa_n(hh)*sum( ht.*exp( -2*pi*1i*(ka_kx(hh)*ra_rx+ka_ky(hh)*ra_ry+ka_kz(hh)*ra_rz)-bf*s2(hh) ) );
end
k = kfactor(Fobs,Fcalc);
Fcalc = Fcalc*k;
dd1 = sum( abs( Fobs(:)-Fcalc(:) ).^2 );

Fsub = (Fcalc-Fobs);
parfor hh = 1:M
Cexp = -ht(hh)*s2.*fa_n.*exp( -2*pi*1i*(ka_kx*ra_rx(hh)+ka_ky*ra_ry(hh)+ka_kz*ra_rz(hh))-bf(hh)*s2 ).*conj(Fsub);
dbf(hh) = 2*real(sum(Cexp));
end

parfor hh = 1:M
Cexp = fa_n.*exp( -2*pi*1i*(ka_kx*ra_rx(hh)+ka_ky*ra_ry(hh)+ka_kz*ra_rz(hh))-bf(hh)*s2 ).*conj(Fsub);
dht(hh) = 2*real(sum(Cexp));
end

for i = 1:Num_type
dd2(i) = sum(dbf(atoms==i));
dd2(Num_type+i) = sum(dht(atoms==i));
end
end