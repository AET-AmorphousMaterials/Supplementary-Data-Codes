function Fcalc = BHfactorFit2(para,xdata)
%fprintf('call BHfactoFit\n');

fa_avg = xdata.fa_avg;
q2 = xdata.q2;
F1 = xdata.F1;
F2 = xdata.F2;
F3 = xdata.F3;

b_vals =  para(1,:);
h_vals =  para(2,:);

g1 =  F1.*exp(-b_vals(1)*q2); 
g2 =  F2.*exp(-b_vals(2)*q2); 
g3 =  F3.*exp(-b_vals(3)*q2); 

Fcalc = fa_avg.*(h_vals(1)*g1 + h_vals(2)*g2 + h_vals(3)*g3);

Fcalc=[real(Fcalc);imag(Fcalc)];
% Fcalc = zeros(L,1);
% for hh = 1:L
%     Fcalc(hh) = fa_avg(hh)*sum( ht.*exp( -2*pi*1i*(kx(hh)*ra.rx+ky(hh)*ra.ry+kz(hh)*ra.rz)-bf*q2(hh) ) );
% end
% R_factor = norm(Fcalc - Fcal,'fro')/norm(Fcalc,'fro')

% k = kfactor(Fobs,Fcalc);
% Fcalc = Fcalc*k;
% Fcalc=Fcalc(:);
end

