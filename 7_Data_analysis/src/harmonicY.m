function Y = harmonicY(n,m,th,phi,varargin)
%HARMONICY  Spherical harmonic function.
%   % th - [0 pi] phi - [0 2*pi]
%   Y = HARMONICY(N,M,TH,PHI) computes the surface spherical harmonic of
%   degree N and order M, evaluated at each element of inclination TH and
%   azimuth PHI. N and M must be scalar integers where M <= abs(N). TH and
%   PHI must be arrays of equal size.
%
%   Y = HARMONICY(N,M,TH,PHI,R) computes the solid spherical harmonic of
%   degree N and order M, evaluated at each element of inclination TH,
%   azimuth PHI and radius R. N and M must be scalar integers where M <=
%   abs(N). TH, PHI and R must be arrays of equal size.
%
%   Y = HARMONICY(__,Name,Value) specifies additional options using one or
%   more name-value pair arguments. Valid options are:
%
%     - 'type' specifies whether to compute the complex spherical harmonics
%       or their real (cosine) or imaginary (sine) parts. Valid values are
%       'complex' (default), 'real' and 'imag'.
%
%     - 'norm' specifies whether the result of the computation is to be
%       normalized. The normalization coefficient is chosen so as to ensure
%       that the spherical harmonics are orthonormal. Valid values are true
%       (default), false.
%       
%     - 'phase' specifies whether to include the Condon-Shortley phase
%       term. This term is not strictly necessary but may be useful in
%       quantum mechanics applications. Valid values are true (default),
%       false.
%
%   See also LEGENDRE.

%   Copyright (C) 2018 Javier Montalt Tordera.

% first argument in varargin to be part of name,value pairs
s = 1;

if nargin > 4
    % if more than four arguments were passed, check varargin{1}
    if ischar(varargin{1})
        % if char array, this is a name and r is empty
        r = [];
    else
        % if not, this should be r
        r = varargin{1};
        % and name,value pairs, if any, start at the second element
        s = 2;
    end
else
    % if only four arguments were passed, r is empty
    r = [];
end

% if r is not empty, we're computing a solid harmonic
issolid = ~isempty(r);

% check inputs
validateattributes(n,{'numeric'},{'scalar','nonempty','integer','nonnegative'},1);
validateattributes(m,{'numeric'},{'scalar','nonempty','integer','>=',-n,'<=',n},2);
validateattributes(th,{'numeric'},{},3);
validateattributes(phi,{'numeric'},{'size',size(th)},4);
if issolid
    validateattributes(r,{'numeric'},{'size',size(th)},5);
end

% parse optional inputs
ip = inputParser;
addParameter(ip,'type','complex',@(x)any(validatestring(x,{'real','imag','complex'})));
addParameter(ip,'norm',true,@(x)islogical(x));
addParameter(ip,'phase',true,@(x)islogical(x));
parse(ip,varargin{s:end});
pp = ip.Results;

% check if m is odd and negative
isoddm = mod(m,2) == 1;
isnegm = m < 0;

% if m is negative, compute the symmetric function where m > 0
if isnegm
    m = abs(m);
end

% to linear arrays
sz = size(th);
th = th(:);
phi = phi(:);
r = r(:);

% normalization factor
if pp.norm
    a = (2*n+1)/(4*pi);
    b = factorial(n-m)/factorial(n+m);
    C = sqrt(a*b);
else
    C = 1;
end

% associated Legendre function
P = legendre(n,cos(th));
P = P(abs(m)+1,:)';

switch pp.type
    case 'real',    E = cos(m*phi);
    case 'imag',    E = sin(m*phi);
    case 'complex', E = exp(1i*m*phi);
end
        
% surface spherical harmonics
Y = C * P .* E;

% solid spherical harmonics
if issolid
    Y = Y .* r.^n;
end

% include Condon-Shortley phase term
% this term is in fact included in MATLAB's associated Legendre function
% definition, so what we are doing here is removing it if not requested
if ~pp.phase && isoddm
    Y = -Y;
end

% if m was negative
if isnegm
    Y = conj(Y);
    if isoddm
        Y = -Y;
    end
end

% reshape to original size
Y = reshape(Y,sz);

end