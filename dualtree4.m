function [A, D] = dualtree4(f, level, useDouble, varargin)
% 4-D Dual-Tree Complex Wavelet Transform
%   [A,D] = DUALTREE4(F) returns the 4-D dual-tree complex wavelet
%   transform of F at the maximum level, floor(log2(min(size(F)))). F is a
%   real-valued 4-D array (X-by-Y-by-Z-by-T) where all dimensions (X,Y,Z,T)
%   must be even and greater than or equal to 4. For now, only the near
%   symmetric biorthogonal wavelet filter with lengths 5 (scaling filter)
%   and 7 (wavelet filter) is used for level 1 and the orthogonal Q-shift
%   Hilbert wavelet filter pair of length 10 is used for levels greater
%   than or equal to 2. A is the matrix of real-valued final-level
%   scaling (lowpass) coefficients. D is a 1-by-L cell array of wavelet
%   coefficients, where L is the level of the transform. There are 120
%   wavelet subbands in the 4-D dual-tree transform at each level. The
%   wavelet coefficients are complex-valued. 
%
%   [A,D] = DUALTREE4(X,LEVEL) obtains the 4-D dual-tree transform down to
%   LEVEL. LEVEL is a positive integer greater than or equal to 2 and less
%   than or equal to floor(log2(min(size(F))).
%
%   USEDOUBLE toggles between double precision arrays (useDouble = 1) and
%   singles precision (useDouble = 0). Double precision is used by default.
%
%   "ExcludeL1" excludes the first level detail coefficients and only the
%   lowpass filter is used. If this option is used a perfect reconstruction
%   is no longer possible but both the forward and inverse transform become
%   computationally more efficient.
%
%   This code is heavily based on the DUALTREE3-function.
%
%   Tommi Heikkilä
%   Created 12.5.2020
%   Last edited 25.1.2021

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

if nargin < 3
    useDouble = 1;
end

% Check for 'ExcludeLeve1Details' or "ExcludeLevel1Details"
validopts = ["ExcludeL1","IncludeL1"];
defaultopt = "IncludeL1";
[opt, varargin] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);

% Case f to double or single
if useDouble
    f = double(f);
else
    f = single(f);
end

% Obtain the first-level analysis filter and q-shift filters
load(char(params.Faf));
load(char(params.af));
if ~useDouble
    LoD = single(LoD);
    HiD = single(HiD);
    LoDa = single(LoDa);
    HiDa = single(HiDa);
end

level = params.level;

% First level filters
h0o = LoD;
h1o = HiD;

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
    A = level1NoHighpass(f,h0o);
    D{1} = [];
else
    [A,D{1}] = level1Highpass(f,h0o,h1o);
end

lev = 2;
% For levels two and up, we use the Qshift filters
while lev <= level    
    [A,D{lev}] = level2Analysis(A,h0a,h1a,h0b,h1b);
    lev = lev+1;
end
end

%------------------------------------------------------------------------
function y = columnFilter(x,h)
% Filter the columns of x with h. This function does not decimate the
% output.

% This determines the symmetric extension of the matrix
L = length(h);
M = fix(L/2);

x = wextend('ar','sym',x,M);
y = conv2(x,h(:),'valid');
end

%------------------------------------------------------------------------
function Z = OddEvenFilter(x,ha,hb)
% ha and hb are identical length (even) filters
[r,c] = size(x); % x is 2-D array
% Even and odd polyphase components of dual-tree filters
haOdd = ha(1:2:end);
haEven = ha(2:2:end);
hbOdd = hb(1:2:end);
hbEven = hb(2:2:end);
r2 = r/2;           % NOTE: THE NUMBER OF ROWS IS HALVED!
Z = zeros(r2,c);
M = length(ha);
% Set up vector for indexing into the matrix
idx = uint8(6:4:r+2*M-2);
matIdx = wextend('ar','sym',(uint8(1:r))',M);

% Now perform the filtering
if dot(ha,hb) > 0
    s1 = uint8(1:2:r2);
    s2 = s1 + 1;
else
    s2 = uint8(1:2:r2);
    s1 = s2 + 1;
end
Z(s1,:) = conv2(x(matIdx(idx-1),:),haOdd(:),'valid') + conv2(x(matIdx(idx-3),:),haEven(:),'valid');
Z(s2,:) = conv2(x(matIdx(idx),:),hbOdd(:),'valid') + conv2(x(matIdx(idx-2),:),hbEven(:),'valid');
end

%-------------------------------------------------------------------------
function A = level1NoHighpass(x,h0o)
% This function is called if the user specified "excludeL1"
sx = size(x);

% Filter dimensions 4 and 3
for rowidx = 1:sx(1)
   for colidx = 1:sx(2)
       % Select one 2D array and permute it so that the 4th dimension
       % forms the columns.
       %                                              t z x y
       y = columnFilter(permute(x(rowidx,colidx,:,:),[4,3,1,2]),h0o);
       % Transpose and filter dimension 3.
       x(rowidx,colidx,:,:) = reshape(columnFilter(y.',h0o),[1,1,sx(3),sx(4)]);
   end
end

% Filter dimensions 2 and 1
for sliceidx = 1:sx(3)
   for timeidx = 1:sx(4)
       % Select one 2D array and transpose it so that the 2nd dimension
       % forms the columns.
       y = columnFilter(x(:,:,sliceidx,timeidx).',h0o);
       % Transpose and filter dimension 1.
       x(:,:,sliceidx,timeidx) = columnFilter(y.',h0o);
   end
end

A = x;
end

%------------------------------------------------------------------------
function [A,D] = level1Highpass(x,h0o,h1o)
% This function computes first level wavelet coefficients

sx = size(x);
if isa(x, 'single')
    xtmp = zeros(2*sx,'single');
else
    xtmp = zeros(2*sx);
end
sxtmp = size(xtmp);
% sr = size(xtmp)/2; -> sr = sx;

% Note this has been extended to be twice the original input size
s1a = uint8(1:sx(1));
s2a = uint8(1:sx(2));
s3a = uint8(1:sx(3));
s4a = uint8(1:sx(4));

s1b = sx(1)+uint8(1:sx(1));
s2b = sx(2)+uint8(1:sx(2));
s3b = sx(3)+uint8(1:sx(3));
s4b = sx(4)+uint8(1:sx(4));

xtmp(s1a,s2a,s3a,s4a) = x; % First orthant of xtmp
clear x

% Filter dimensions 4 and 3
% Note that we only loop through normal number of rows and columns!
for rowidx = 1:sx(1)
   for colidx = 1:sx(2)
       % Select one 2D array and permute it so that the 4th dimension
       % forms the columns.
       %                                        t z x y
       y = permute(xtmp(rowidx,colidx,s3a,s4a),[4,3,1,2]); 
       %       Lowpass              Highpass
       y = [columnFilter(y,h0o); columnFilter(y,h1o)].'; % Combine and transpose
       sy = cat(2,[1,1],size(y)); % Consider size as a part of 4-D array
       
       % 3rd dimension
       xtmp(rowidx,colidx,s3a,:) = reshape(columnFilter(y,h0o),sy); % Lowpass
       xtmp(rowidx,colidx,s3b,:) = reshape(columnFilter(y,h1o),sy); % Highpass
   end
end

% Filter dimensions 2 and 1
% Note that we loop through doubled number of slices and time steps!
for sliceidx = 1:sxtmp(3)
    for timeidx = 1:sxtmp(4)
        % Select one 2D array and permute it so that the 2nd dimension
        % forms the columns.
        %                                           y x z t
        y = permute(xtmp(s1a,s2a,sliceidx,timeidx),[2,1,3,4]); 
        %       Lowpass              Highpass
        y = [columnFilter(y,h0o); columnFilter(y,h1o)].'; % Combine and transpose
        sy = cat(2,[1,1],size(y)); % Consider size as a part of 4-D array
        
        % 1st dimension
        xtmp(s1a,:,sliceidx,timeidx) = reshape(columnFilter(y,h0o),sy); % Lowpass
        xtmp(s1b,:,sliceidx,timeidx) = reshape(columnFilter(y,h1o),sy); % Highpass
    end
end

% Note in listing the subbands the order is reversed compared to what was
% done previously, i.e. 1st dimension, then 2nd, then 3rd and then 4th.
A = xtmp(s1a,s2a,s3a,s4a);                                  % LLLL
% Form the eight complex wavelets for 4^2-1 = 15 subbands for a total of
% 15*8 = 120 sets of coefficients per level.
D = cat(5, cube2complex(xtmp(s1b,s2a,s3a,s4a)),...          % HLLL
            cube2complex(xtmp(s1a,s2b,s3a,s4a)),...         % LHLL
            cube2complex(xtmp(s1b,s2b,s3a,s4a)),...         % HHLL
            cube2complex(xtmp(s1a,s2a,s3b,s4a)),...         % LLHL
            cube2complex(xtmp(s1b,s2a,s3b,s4a)),...         % HLHL
            cube2complex(xtmp(s1a,s2b,s3b,s4a)),...         % LHHL
            cube2complex(xtmp(s1b,s2b,s3b,s4a)),...         % HHHL
            cube2complex(xtmp(s1a,s2a,s3a,s4b)),...         % LLLH
            cube2complex(xtmp(s1b,s2a,s3a,s4b)),...         % HLLH
            cube2complex(xtmp(s1a,s2b,s3a,s4b)),...         % LHLH
            cube2complex(xtmp(s1b,s2b,s3a,s4b)),...         % HHLH
            cube2complex(xtmp(s1a,s2a,s3b,s4b)),...         % LLHH
            cube2complex(xtmp(s1b,s2a,s3b,s4b)),...         % HLHH
            cube2complex(xtmp(s1a,s2b,s3b,s4b)),...         % LHHH
            cube2complex(xtmp(s1b,s2b,s3b,s4b)));           % HHHH
end

%-------------------------------------------------------------------
function z = cube2complex(x)
% Form the complex-valued subbands
J = 1/2*[1 1j];
% Form 16 building blocks P(si) from a tree-like structure. Even indexing
% corresponds real valued part (tree a), odd to imaginary part (tree b).
%
% Even number of imaginary or real parts corresponds to real.
% Odd number of imaginary or real parts corresponds to imaginary.
%            X       Y       Z       T
Paaaa = x(2:2:end,2:2:end,2:2:end,2:2:end); % Re
Paaab = x(2:2:end,2:2:end,2:2:end,1:2:end); % Im
Paaba = x(2:2:end,2:2:end,1:2:end,2:2:end); % Im
Paabb = x(2:2:end,2:2:end,1:2:end,1:2:end); % Re
Pabaa = x(2:2:end,1:2:end,2:2:end,2:2:end); % Im
Pabab = x(2:2:end,1:2:end,2:2:end,1:2:end); % Re
Pabba = x(2:2:end,1:2:end,1:2:end,2:2:end); % Re
Pabbb = x(2:2:end,1:2:end,1:2:end,1:2:end); % Im
Pbaaa = x(1:2:end,2:2:end,2:2:end,2:2:end); % Im
Pbaab = x(1:2:end,2:2:end,2:2:end,1:2:end); % Re
Pbaba = x(1:2:end,2:2:end,1:2:end,2:2:end); % Re
Pbabb = x(1:2:end,2:2:end,1:2:end,1:2:end); % Im
Pbbaa = x(1:2:end,1:2:end,2:2:end,2:2:end); % Re
Pbbab = x(1:2:end,1:2:end,2:2:end,1:2:end); % Im
Pbbba = x(1:2:end,1:2:end,1:2:end,2:2:end); % Im
Pbbbb = x(1:2:end,1:2:end,1:2:end,1:2:end); % Re
clear x

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
z = cat(5,O1,O2,O3,O4,O5,O6,O7,O8);
end

%-------------------------------------------------------------------------
function [A,D] = level2Analysis(x,h0a,h1a,h0b,h1b)
% This the analysis bank for levels >= 2, here we require the four qshift
% filters
% First we want to guarantee that the input LLLL image is divisible by
% four in each dimension

LLLLsize = size(x);
if any(rem(LLLLsize,4))
    x = paddata(x);
    % Now get size of extended x
    LLLLsize = size(x);
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

%
% Filter dimensions 4 and 3
% Note that we loop through full number of rows and columns!
for rowidx = 1:LLLLsize(1)
   for colidx = 1:LLLLsize(2)
       % Select one 2D array and permute it so that the 4th dimension
       % forms the columns.
       %                                 t z x y
       y = permute(x(rowidx,colidx,:,:),[4,3,1,2]); 
       %       Lowpass                   Highpass
       y = [OddEvenFilter(y,h0b,h0a); OddEvenFilter(y,h1b,h1a)].'; % Combine and transpose
       sy = [1,1,size(y,1)/2,size(y,2)]; % Consider size as a part of 4-D array
       
       % 3rd dimension
       x(rowidx,colidx,s3a,:) = reshape(OddEvenFilter(y,h0b,h0a),sy); % Lowpass
       x(rowidx,colidx,s3b,:) = reshape(OddEvenFilter(y,h1b,h1a),sy); % Highpass
   end
end

% Filter dimensions 2 and 1
% Note that we loop through full number of slices and time steps!
for sliceidx = 1:LLLLsize(3)
    for timeidx = 1:LLLLsize(4)
        % Select one 2D array and permute it so that the 2nd dimension
        % forms the columns.
        %                                    y x z t
        y = permute(x(:,:,sliceidx,timeidx),[2,1,3,4]); 
        %       Lowpass                   Highpass
        y = [OddEvenFilter(y,h0b,h0a); OddEvenFilter(y,h1b,h1a)].'; % Combine and transpose
        sy = [1,1,size(y,1)/2,size(y,2)]; % Consider size as a part of 4-D array
        
        % 1st dimension
        x(s1a,:,sliceidx,timeidx) = reshape(OddEvenFilter(y,h0b,h0a),sy); % Lowpass
        x(s1b,:,sliceidx,timeidx) = reshape(OddEvenFilter(y,h1b,h1a),sy); % Highpass
    end
end


% Form the eight complex wavelets for 4^2-1 = 15 subbands for a total of
% 15*8 = 120 sets of coefficients per level.
D = cat(5, cube2complex(x(s1b,s2a,s3a,s4a)),...             % HLLL
            cube2complex(x(s1a,s2b,s3a,s4a)),...            % LHLL
            cube2complex(x(s1b,s2b,s3a,s4a)),...            % HHLL
            cube2complex(x(s1a,s2a,s3b,s4a)),...            % LLHL
            cube2complex(x(s1b,s2a,s3b,s4a)),...            % HLHL
            cube2complex(x(s1a,s2b,s3b,s4a)),...            % LHHL
            cube2complex(x(s1b,s2b,s3b,s4a)),...            % HHHL
            cube2complex(x(s1a,s2a,s3a,s4b)),...            % LLLH
            cube2complex(x(s1b,s2a,s3a,s4b)),...            % HLLH
            cube2complex(x(s1a,s2b,s3a,s4b)),...            % LHLH
            cube2complex(x(s1b,s2b,s3a,s4b)),...            % HHLH
            cube2complex(x(s1a,s2a,s3b,s4b)),...            % LLHH
            cube2complex(x(s1b,s2a,s3b,s4b)),...            % HLHH
            cube2complex(x(s1a,s2b,s3b,s4b)),...            % LHHH
            cube2complex(x(s1b,s2b,s3b,s4b)));              % HHHH

% This subband returned as a matrix because only the coarsest
% resolution is retained.
A = x(s1a,s2a,s3a,s4a);                                     % LLLL
end

%-------------------------------------------------------------------------
function x = paddata(x)
% Pad data if necessary
sx = size(x);
if rem(sx(1),4)
    x = cat(1,x(1,:,:,:),x,x(end,:,:,:));
end
if rem(sx(2),4)
    x = cat(2,x(:,1,:,:),x,x(:,end,:,:));
end
if rem(sx(3),4)
    x = cat(3,x(:,:,1,:),x,x(:,:,end,:));
end
if rem(sx(4),4)
    x = cat(4,x(:,:,:,1),x,x(:,:,:,end));
end
end
