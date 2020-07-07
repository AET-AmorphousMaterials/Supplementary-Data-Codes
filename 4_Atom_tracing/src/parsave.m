function parsave(varargin)
savefile = varargin{1}; % first input argument
for i = 2:2:nargin
    savevar.(varargin{i}) = varargin{i+1}; % other input arguments
end
save(savefile,'-struct','savevar')
end