function fixedfa =  make_fixedfa_man3d(sizeXY, Res, Z_arr)
    finalvol_summed = zeros(sizeXY);

    kx = 1:size(finalvol_summed,1);
    ky = 1:size(finalvol_summed,2);
    kz = 1:size(finalvol_summed,3);

    MultF_X = 1/(length(kx)*Res);
    MultF_Y = 1/(length(ky)*Res);
    MultF_Z = 1/(length(kz)*Res);

    CentPos = round((size(finalvol_summed)+1)/2);
    [KX, KY, KZ] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y, (kz-CentPos(3))*MultF_Z);
    q2 = KX.^2 + KY.^2 +KZ.^2;
    clear KX KY KZ

    fixedfa_arr = zeros(numel(Z_arr),numel(q2));
    for i = 1:numel(Z_arr)
        fixedfa_arr(i,:) = fatom_vector(sqrt(q2),Z_arr(i));
    end
    fixedfa = mean(fixedfa_arr,1);
    fixedfa = reshape(fixedfa, sizeXY);
end