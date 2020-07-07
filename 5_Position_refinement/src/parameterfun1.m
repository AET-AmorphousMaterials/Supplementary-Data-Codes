function [dd1, dd2, k] = parameterfun1(x,ht,bf,ra,fa,ka,L,Fobs,Fcalc,atoms)

bf(atoms==1) = x(1); bf(atoms==2) = x(2);
ht(atoms==1) = x(3); ht(atoms==2) = x(4);

M = size(ra.rx,2);
dbf = zeros(1,M);
dht = zeros(1,M);
s2 = ka.kx.^2+ka.ky.^2+ka.kz.^2;
fa_3 = fa(3,:);

parfor hh = 1:L
Fcalc(hh) = fa_3(hh)*sum( ht.*exp( -2*pi*1i*(ka.kx(hh)*ra.rx+ka.ky(hh)*ra.ry+ka.kz(hh)*ra.rz)-bf*s2(hh) ) );
end
k = kfactor(Fobs,Fcalc);
Fcalc = Fcalc*k;
dd1 = sum( abs( Fobs(:)-Fcalc(:) ).^2 );



Fsub = (Fcalc-Fobs);
parfor hh = 1:M
Cexp = -ht(hh)*s2.*fa_3.*exp( -2*pi*1i*(ka.kx*ra.rx(hh)+ka.ky*ra.ry(hh)+ka.kz*ra.rz(hh))-bf(hh)*s2 ).*conj(Fsub);
dbf(hh) = 2*real(sum(Cexp));
end
dd2(1) = sum(dbf(atoms==1));
dd2(2) = sum(dbf(atoms==2));



parfor hh = 1:M
Cexp = fa_3.*exp( -2*pi*1i*(ka.kx*ra.rx(hh)+ka.ky*ra.ry(hh)+ka.kz*ra.rz(hh))-bf(hh)*s2 ).*conj(Fsub);
dht(hh) = 2*real(sum(Cexp));
end
dd2(3) = sum(dht(atoms==1));
dd2(4) = sum(dht(atoms==2));