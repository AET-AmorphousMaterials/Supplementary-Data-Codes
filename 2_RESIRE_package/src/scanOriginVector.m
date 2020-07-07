function vectorsOrigin = scanOriginVector(zeroVector,AnglesZYX)

Z = [0 0 1]';
Y = [0 1 0]';
X = [1 0 0]';

% MAT1 = MatrixQuaternionRot([0 0 1],ROT(1));
% MAT2 = MatrixQuaternionRot([0 1 0],ROT(2));
% RotVeca = MAT2*MAT1*(-X);
%  
% MAT1 = MatrixQuaternionRot([0 0 1],ROT(3));
% MAT2 = MatrixQuaternionRot([1 0 0],ROT(4));
% RotVecb = MAT2*MAT1*Y;

% AnglesZYX_GENFIRE_noimagerot = zeros(size(AnglesCalib));

if size(zeroVector,1)==1
    zeroVector = zeroVector';
end

vectorsOrigin = zeros(size(AnglesZYX,1),3);
for pj = 1:size(AnglesZYX,1)
    R1 = MatrixQuaternionRot(Z, AnglesZYX(pj,1));
    R2 = MatrixQuaternionRot(Y, AnglesZYX(pj,2));
    R3 = MatrixQuaternionRot(X, AnglesZYX(pj,3));
    ROTMAT = (R1*R2*R3); %applied to coordinates
    vectorsOrigin(pj,:)=(ROTMAT'*zeroVector)'; % applied to object
    
%     ROTMnew = ROTMAT*RotZero;
%     AnglesZYX_GENFIRE_noimagerot(pj,:) = get_ZYX_Euler_angles_from_matrix(ROTMnew');
end
