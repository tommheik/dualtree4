function [A, D] = dualtree4(f, level, varargin)
% 4-D Dual-Tree Complex Wavelet Transform
%   [A,D] = DUALTREE4(F) returns the 4-D dual-tree complex wavelet
%   transform of F at the maximum level, floor(log2(min(size(F)))). F is a
%   real-valued 4-D array (X-by-Y-by-Z-by-T) where all dimensions (X,Y,Z,T)
%   must be even and greater than or equal to 4. For now, only the near
%   symmetric biorthogonal wavelet filter with lengths 5 (scaling filter)
%   and 7 (wavelet filter) is used for level 1 and the orthogonal Q-shift
%   Hilbert wavelet filter pair of length 10 is used for levels greater
%   than or equal to 2.
%
%   A is the matrix of complex or real-valued final-level scaling (lowpass)
%   coefficients. By default complex values are used.
%
%   D is a 1-by-L cell array of wavelet coefficients, where L is the level
%   of the transform. There are 120 wavelet subbands in the 4-D dual-tree
%   transform at each level. The wavelet coefficients are complex-valued. 
%
%   [A,D] = DUALTREE4(X,LEVEL) obtains the 4-D dual-tree transform down to
%   LEVEL. LEVEL is a positive integer greater than or equal to 2 and less
%   than or equal to floor(log2(min(size(F))).
%
%   "ExcludeL1" excludes the first level detail coefficients and only the
%   lowpass filter is used. If this option is used a perfect reconstruction
%   is no longer possible but both the forward and inverse transform become
%   computationally more efficient.
%
%   This code is heavily based on the DUALTREE3-function.
%
%   Tommi Heikkilä
%   University of Helsinki, Dept. of Mathematics and Statistics
%   Created 12.5.2020
%   Last edited 7.5.2021

% Ensure the input is numeric, real, and that it is four-dimensional
validateattributes(f,{'numeric'},{'real','nonempty','finite','ndims',4},...
    'DUALTREE4','F');

% Check that all the dimensions of x are even and every dimension is
% greater than or equal to 4
origsizedata = size(f);
if any(rem(origsizedata,2)) || any(origsizedata < 4)
    error('Object dimensions are incompatible');
end
% Set up defaults for the first-level biorthogonal filter and subsequent
% level orthogonal Hilbert filters
% NOTE: CUSTOM FILTER LENGTHS NOT YET IMPLEMENTED!
maxlev = floor(log2(min(size(f))));
params.Faf = "nearsym5_7";
params.af = 10;
params.level = maxlev;
params.af = deblank(['qshift' num2str(params.af)]);
if nargin >= 2
validateattributes(level,{'numeric'},...
        {'integer','scalar','<=',maxlev,'>=',1},'DUALTREE4','LEVEL');
    params.level = level;
end

% Use double precision if f is double, otherwise use single
useDouble = isa(f,'double');

% Check for 'ExcludeLeve1Details' or "ExcludeLevel1Details"
validopts = ["ExcludeL1","IncludeL1"];
defaultopt = "IncludeL1";
[opt, varargin] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);

% Check for 'realA' or "realA", i.e. whether the approximation
% coefficients are returns as complex coefficients.
validopts = ["realA","complexA"];
defaultopt = "complexA";
[Atype, ~] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);

% Case f to double or single
if ~useDouble
    f = single(f);
end

% Obtain the first-level analysis filter and q-shift filters
load(char(params.Faf),'LoD','HiD');
load(char(params.af),'LoD*','HiD*');
if ~useDouble
    LoD = single(LoD);
    HiD = single(HiD);
    LoDa = single(LoDa);
    HiDa = single(HiDa);
end

level = params.level;

% First level filters
h0 = LoD;
h1 = HiD;

% Filters for levels >= 2
h0a = LoDa;
h1a = HiDa;

% Normalize analysis filters
hscale = 1 / norm(h0a,2);
% Tree A analysis filters
h0a = h0a.* hscale;
h1a = h1a.* hscale;

% Tree B analysis filters
h0b = h0a(end:-1:1);
h1b = h1a(end:-1:1);

% Debug cleanup
clear Hi* Lo*

% Allocate array for wavelet coefficients
D = cell(level,1);

% Level 1 filtering. We can omit the highest level
% details if needed
if strcmpi(opt,"ExcludeL1")
    A = level1NoHighpass(f,h0);
    D{1} = [];
else
    [A,D{1}] = level1Highpass(f,h0,h1);
end

lev = 2;
% For levels two and up, we use the Qshift filters
while lev <= level    
    [A,D{lev}] = level2Analysis(A,h0a,h1a,h0b,h1b);
    lev = lev+1;
end
if strcmpi(Atype,"complexA") % Check if A is returned complex valued
    A = cube2complex(A);
end
end

%------------------------------------------------------------------------
function X = oneDFilter(X,h,dim)
% Filter one dimension of X with h, where h is a column vector. The output
% is NOT downsampled

% Determine symmetric extension amount
h = h(:);
lh = length(h);
a = fix(lh/2);

% Permute h so that it is "dim-dimensional" vector
switch dim
    case 1
        % Do nothing
    case 2
        h = h';
    case 3
        h = reshape(h,1,1,[]);
    case 4
        h = reshape(h,1,1,1,[]);
end

% Extend X and convolve
Y = wextend4D(X,a,dim);
X = convn(Y,h,'valid');

end

%------------------------------------------------------------------------
function Y = wextend4D(X,a,dim)
% 4D version of wextend using symmetric half-point extension on rows
lx = size(X,dim);
% We get the indicies from the original wextend
i = wextend('ac','sym',1:lx,a);
% Create cell array where i is placed on dim'th place
I = cell(1,4); I(:) = {':'};
I{dim} = i;
% Extend dim'th direction using i
Y = X(I{:});
end

%------------------------------------------------------------------------
function Z = OddEvenFilter4D(X,ha,hb,dim)
% Dual filter scheme where the convolutions using filters ha and hb are
% interlaced.
% THIS FUNCTION HALVES THE SIZE OF THE CONVOLVED DIRECTION!
% This is because the input for level 2 and up has NOT been downsampled
% yet and hence it is twice the size it should be.

% ha and hb are identical length (even) filters
% Case to column vectors
ha = ha(:);
hb = hb(:);
M = length(ha);

% permute filters
switch dim
    case 1
        % Do nothing
    case 2
        ha = ha';
        hb = hb';
    case 3
        ha = reshape(ha,1,1,[]);
        hb = reshape(hb,1,1,[]);
    case 4
        ha = reshape(ha,1,1,1,[]);
        hb = reshape(hb,1,1,1,[]);
end

szX = size(X); % X is 4-D array
% Even and odd polyphase components of dual-tree filters
haOdd = ha(1:2:end);
haEven = ha(2:2:end);
hbOdd = hb(1:2:end);
hbEven = hb(2:2:end);

od = szX(dim)/2;           % NOTE: operated dimension is halved
szZ = szX; szZ(dim) = od;

% Initialize
Z = zeros(szZ,class(X));
% Set up vector for indexing into the matrix
skipInd = uint8(6:4:szX(dim)+2*M-2);
extIdx = wextend('ac','sym',(uint8(1:szX(dim))),M);

% Now perform the filtering
if dot(ha,hb) > 0
    s1 = uint8(1:2:od); % Odd values
    s2 = s1 + 1; % Even values
else
    s2 = uint8(1:2:od); % Odd values
    s1 = s2 + 1; % Even values
end
J = cell(1,4); J(:) = {':'};
% Filter with hb
% Create cell array where correct indices are placed on dim'th place
Iodd = J;
Ieven = J;
Iodd{dim}   = extIdx(skipInd-1);
Ieven{dim}  = extIdx(skipInd-3);
J{dim}      = s1;

Z(J{:}) = convn(X(Iodd{:}),hbOdd,'valid') ...
    + convn(X(Ieven{:}),hbEven,'valid');

% Filter with ha
% Change the cell array accordingly
Iodd{dim}   = extIdx(skipInd);
Ieven{dim}  = extIdx(skipInd-2);
J{dim}      = s2;

Z(J{:}) = convn(X(Iodd{:}),haOdd,'valid') ...
    + convn(X(Ieven{:}),haEven,'valid');
end

%-------------------------------------------------------------------------
function A = level1NoHighpass(x,h0)
% This function is called if the user specified "excludeL1"

% Filter dimension 4
y = oneDFilter(x,h0,4);

% Filter dimension 3
x = oneDFilter(y,h0,3);

% Filter dimension 2
y = oneDFilter(x,h0,2);

% Filter dimension 1
A = oneDFilter(y,h0,1);
end

%-------------------------------------------------------------------------
function [A,D] = level1Highpass(X,h0,h1)
% This function computes first level wavelet coefficients

sX = size(X);

% Note this has been extended to be twice the original input size
s2a = uint8(1:sX(2));
s3a = uint8(1:sX(3));
s4a = uint8(1:sX(4));

s2b = sX(2)+uint8(1:sX(2));
s3b = sX(3)+uint8(1:sX(3));
s4b = sX(4)+uint8(1:sX(4));

% It is faster to work with two smaller arrays than one big array

% Filter dimension 4
Yl = oneDFilter(X,h0,4); % Lowpass
Yh = oneDFilter(X,h1,4); % Highpass
Y = cat(4,Yl,Yh);

% Filter dimension 3
Xl = oneDFilter(Y,h0,3); % Lowpass
Xh = oneDFilter(Y,h1,3); % Highpass
X = cat(3,Xl,Xh);

% Filter dimension 2
Yl = oneDFilter(X,h0,2); % Lowpass
Yh = oneDFilter(X,h1,2); % Highpass
Y = [Yl Yh];

% Filter dimension 1
Xl = oneDFilter(Y,h0,1); % Lowpass
Xh = oneDFilter(Y,h1,1); % Highpass

% Note in listing the subbands the order is reversed compared to what was
% done previously, i.e. 1st dimension, then 2nd, then 3rd and then 4th.
A = Xl(:,s2a,s3a,s4a);                  % LLLL
% Form the eight complex wavelets for 4^2-1 = 15 subbands for a total of
% 15*8 = 120 sets of coefficients per level.
Y1  = cube2complex(Xh(:,s2a,s3a,s4a));  % HLLL
Y2  = cube2complex(Xl(:,s2b,s3a,s4a));  % LHLL
Y3  = cube2complex(Xh(:,s2b,s3a,s4a));  % HHLL
Y4  = cube2complex(Xl(:,s2a,s3b,s4a));  % LLHL
Y5  = cube2complex(Xh(:,s2a,s3b,s4a));  % HLHL
Y6  = cube2complex(Xl(:,s2b,s3b,s4a));  % LHHL
Y7  = cube2complex(Xh(:,s2b,s3b,s4a));  % HHHL
Y8  = cube2complex(Xl(:,s2a,s3a,s4b));  % LLLH
Y9  = cube2complex(Xh(:,s2a,s3a,s4b));  % HLLH
Y10 = cube2complex(Xl(:,s2b,s3a,s4b));  % LHLH
Y11 = cube2complex(Xh(:,s2b,s3a,s4b));  % HHLH
Y12 = cube2complex(Xl(:,s2a,s3b,s4b));  % LLHH
Y13 = cube2complex(Xh(:,s2a,s3b,s4b));  % HLHH
Y14 = cube2complex(Xl(:,s2b,s3b,s4b));  % LHHH
Y15 = cube2complex(Xh(:,s2b,s3b,s4b));  % HHHH

D = cat(5,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y15);          
end

%-------------------------------------------------------------------
function Z = cube2complex(X)
% Form the complex-valued subbands
J = 1/2*[1 1j];
% Form 16 building blocks P(si) from a tree-like structure. Even indexing
% corresponds real valued part (tree a), odd to imaginary part (tree b).
%
% Even number of imaginary or real parts corresponds to real.
% Odd number of imaginary or real parts corresponds to imaginary.
%            X       Y       Z       T
Paaaa = X(2:2:end,2:2:end,2:2:end,2:2:end); % Re
Paaab = X(2:2:end,2:2:end,2:2:end,1:2:end); % Im
Paaba = X(2:2:end,2:2:end,1:2:end,2:2:end); % Im
Paabb = X(2:2:end,2:2:end,1:2:end,1:2:end); % Re
Pabaa = X(2:2:end,1:2:end,2:2:end,2:2:end); % Im
Pabab = X(2:2:end,1:2:end,2:2:end,1:2:end); % Re
Pabba = X(2:2:end,1:2:end,1:2:end,2:2:end); % Re
Pabbb = X(2:2:end,1:2:end,1:2:end,1:2:end); % Im
Pbaaa = X(1:2:end,2:2:end,2:2:end,2:2:end); % Im
Pbaab = X(1:2:end,2:2:end,2:2:end,1:2:end); % Re
Pbaba = X(1:2:end,2:2:end,1:2:end,2:2:end); % Re
Pbabb = X(1:2:end,2:2:end,1:2:end,1:2:end); % Im
Pbbaa = X(1:2:end,1:2:end,2:2:end,2:2:end); % Re
Pbbab = X(1:2:end,1:2:end,2:2:end,1:2:end); % Im
Pbbba = X(1:2:end,1:2:end,1:2:end,2:2:end); % Im
Pbbbb = X(1:2:end,1:2:end,1:2:end,1:2:end); % Re
clear X

% Directionality is obtained by considering the complex conjugates of some
% blocks P. This changes the sign of the imaginary units for different
% directions. For each direction or orthant the real and imaginary blocks
% are combined into a single complex wavelet.
%
% 1st orthant: j1 = j2 = j3 = j4 = +j
O1 = J(1)*(Paaaa-Paabb-Pabab-Pabba-Pbaab-Pbaba-Pbbaa+Pbbbb) + ...
     J(2)*(Paaab+Paaba+Pabaa-Pabbb+Pbaaa-Pbabb-Pbbab-Pbbba);
% 2nd orthant: j1 = -j; j2 = j3 = j4 = +j
O2 = J(1)*(Paaaa-Paabb-Pabab-Pabba+Pbaab+Pbaba+Pbbaa-Pbbbb) + ...
     J(2)*(Paaab+Paaba+Pabaa-Pabbb-Pbaaa+Pbabb+Pbbab+Pbbba);
% 3rd orthant: j2 = -j; j1 = j3 = j4 = +j
O3 = J(1)*(Paaaa-Paabb+Pabab+Pabba-Pbaab-Pbaba+Pbbaa-Pbbbb) + ...
     J(2)*(Paaab+Paaba-Pabaa+Pabbb+Pbaaa-Pbabb+Pbbab+Pbbba);
% 4th orthant: j1 = j2 = -j; j3 = j4 = +j
O4 = J(1)*(Paaaa-Paabb+Pabab+Pabba+Pbaab+Pbaba-Pbbaa+Pbbbb) + ...
     J(2)*(Paaab+Paaba-Pabaa+Pabbb-Pbaaa+Pbabb-Pbbab-Pbbba);
% 5th orthant: j1 = j2 = j4 = +j; j3 = -j
O5 = J(1)*(Paaaa+Paabb-Pabab+Pabba-Pbaab+Pbaba-Pbbaa-Pbbbb) + ...
     J(2)*(Paaab-Paaba+Pabaa+Pabbb+Pbaaa+Pbabb-Pbbab+Pbbba);
% 6th orthant: j1 = j3 = -j; j2 = j4 = +j
O6 = J(1)*(Paaaa+Paabb-Pabab+Pabba+Pbaab-Pbaba+Pbbaa+Pbbbb) + ...
     J(2)*(Paaab-Paaba+Pabaa+Pabbb-Pbaaa-Pbabb+Pbbab-Pbbba);
% 7th orthant: j1 = j4 = +j; j2 = j3 = -j
O7 = J(1)*(Paaaa+Paabb+Pabab-Pabba-Pbaab+Pbaba+Pbbaa+Pbbbb) + ...
     J(2)*(Paaab-Paaba-Pabaa-Pabbb+Pbaaa+Pbabb+Pbbab-Pbbba);
% 8th orthant: j1 = j2 = j3 = -j; j4 = +j
O8 = J(1)*(Paaaa+Paabb+Pabab-Pabba+Pbaab-Pbaba-Pbbaa-Pbbbb) + ...
     J(2)*(Paaab-Paaba-Pabaa-Pabbb-Pbaaa-Pbabb-Pbbab+Pbbba);

% Return all (eight) 4-D objects in one 5-D array.
Z = cat(5,O1,O2,O3,O4,O5,O6,O7,O8);
end

%-------------------------------------------------------------------------
function [A,D] = level2Analysis(X,h0a,h1a,h0b,h1b)
% This the analysis bank for levels >= 2, here we require the four qshift
% filters
% First we want to guarantee that the input LLLL image is divisible by
% four in each dimension

LLLLsize = size(X);
if any(rem(LLLLsize,4))
    X = paddata(X);
    % Now get size of extended X
    LLLLsize = size(X);
end

% These will be integers
sr = LLLLsize/2;
% Set up index vectors for filtering
s1a = uint8(1:sr(1));
s2a = uint8(1:sr(2));
s3a = uint8(1:sr(3));
s4a = uint8(1:sr(4));
s1b = s1a+sr(1);
s2b = s2a+sr(2);
s3b = s3a+sr(3);
s4b = s4a+sr(4);

% We need to keep the input unchanged until both lowpass and highpass
% filters have been used.
Y = zeros(LLLLsize,class(X));

% Filter dimension 4
perm = [4 2 3 1];
Y(:,:,:,s4b) = OddEvenFilter4D(X,h1a,h1b,4); % Highpass
Y(:,:,:,s4a) = OddEvenFilter4D(X,h0a,h0b,4); % Lowpass

% Filter dimension 3
perm = [3 2 1 4];
X(:,:,s3b,:) = OddEvenFilter4D(Y,h1a,h1b,3); % Highpass
X(:,:,s3a,:) = OddEvenFilter4D(Y,h0a,h0b,3); % Lowpass

% Filter dimension 2
perm = [2 1 3 4];
Y(:,s2b,:,:) = OddEvenFilter4D(X,h1a,h1b,2); % Highpass
Y(:,s2a,:,:) = OddEvenFilter4D(X,h0a,h0b,2); % Lowpass

% Filter dimension 1
perm = []; % Same as [1 2 3 4]
X(s1b,:,:,:) = OddEvenFilter4D(Y,h1a,h1b,1); % Highpass
X(s1a,:,:,:) = OddEvenFilter4D(Y,h0a,h0b,1); % Lowpass

% Form the eight complex wavelets for 4^2-1 = 15 subbands for a total of
% 15*8 = 120 sets of coefficients per level.
D = cat(5, cube2complex(X(s1b,s2a,s3a,s4a)),...             % HLLL
            cube2complex(X(s1a,s2b,s3a,s4a)),...            % LHLL
            cube2complex(X(s1b,s2b,s3a,s4a)),...            % HHLL
            cube2complex(X(s1a,s2a,s3b,s4a)),...            % LLHL
            cube2complex(X(s1b,s2a,s3b,s4a)),...            % HLHL
            cube2complex(X(s1a,s2b,s3b,s4a)),...            % LHHL
            cube2complex(X(s1b,s2b,s3b,s4a)),...            % HHHL
            cube2complex(X(s1a,s2a,s3a,s4b)),...            % LLLH
            cube2complex(X(s1b,s2a,s3a,s4b)),...            % HLLH
            cube2complex(X(s1a,s2b,s3a,s4b)),...            % LHLH
            cube2complex(X(s1b,s2b,s3a,s4b)),...            % HHLH
            cube2complex(X(s1a,s2a,s3b,s4b)),...            % LLHH
            cube2complex(X(s1b,s2a,s3b,s4b)),...            % HLHH
            cube2complex(X(s1a,s2b,s3b,s4b)),...            % LHHH
            cube2complex(X(s1b,s2b,s3b,s4b)));              % HHHH

% This subband returned as a matrix because only the coarsest
% resolution is retained.
A = X(s1a,s2a,s3a,s4a);                                     % LLLL
end

%-------------------------------------------------------------------------
function X = paddata(X)
% Pad data if necessary
sx = size(X);
if rem(sx(1),4)
    X = cat(1,X(1,:,:,:),X,X(end,:,:,:));
end
if rem(sx(2),4)
    X = cat(2,X(:,1,:,:),X,X(:,end,:,:));
end
if rem(sx(3),4)
    X = cat(3,X(:,:,1,:),X,X(:,:,end,:));
end
if rem(sx(4),4)
    X = cat(4,X(:,:,:,1),X,X(:,:,:,end));
end
end
