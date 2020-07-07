%% RESIRE reconstruction
function obj = reconstruct_control(obj)
if strcmp(obj.method,'FST')
    obj = runGridding(obj); 
    obj = reconstruct(obj);
elseif strcmp(obj.method,'Radon')
    if sum(sum(abs(obj.InputAngles(:,[1 3]))))==0 && sum(sum(abs(obj.InputAngles(:,2))))~=0
        obj = reconstruct3D2D(obj);
    else
        error('RESIRE: Only single y-tilt implemented!')
    end
else
    error('RESIRE: Unrecognized forward projection method!')
end
end