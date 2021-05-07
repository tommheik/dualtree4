function rec = idualtree4(A, D, varargin)
% 4-D Dual-Tree Wavelet Reconstruction
%   REC = IDUALTREE4(A,D) returns the inverse 4-D dual-tree complex
%   wavelet transform of the final-level approximation coefficients, A, and
%   cell array of wavelet coefficients, D. A and D are outputs of
%   DUALTREE4.
%
%   REC = IDUALTREE4(A,D,orgSize) computes the final synthesis so that the
%   size of REC matches the size of the original object (orgSize). This is
%   sometimes needed if the 1st level detail coefficients were excluded
%   during analysis ([A,D] = DUALTREE4(...,"ExcludeL1")).
%
%   ADJ = IDUALTREE4(A,D,"adjoint") approximates the adjoint operator of the
%   forward transform. For adjoint the 1st level decomposition filters need 
%   to be changed and a normalization identical to the analysis (DUALTREE4)
%   needs to be used.
%
%   This code is heavily based on the IDUALTREE3-function.
%
%   Tommi Heikkilä
%   University of Helsinki, Dept. of Mathematics and Statistics
%   Created 15.5.2020
%   Last edited 7.5.2021

% Check whether approximation coefficients are stored in real or complex 
% valued format.
realA = isreal(A);
if realA
    % A should be a 4-D real-valued matrix output from dualtree4
    validateattributes(A,{'real','numeric'},{'nonempty','finite','ndims',4},...
    'IDUALTREE4','A');
else
    % A should be a 5-D complex-valued matrix output from dualtree4
    validateattributes(A,{'complex','numeric'},{'nonempty','finite','ndims',5},...
    'IDUALTREE4','A');
end
% D should be a cell array
validateattributes(D,{'cell'},{'nonempty'},'IDUALTREE4','D');

% Obtain the level of the transform
level = length(D);

% Use double precision if needed
useDouble = isa(A, 'double');

% Check for 'adjoint' or "adjoint"
validopts = ["adjoint","inverse"];
defaultopt = "inverse";
global useAdj % Simpler to pass to other functions this way
[useAdj, varargin] = ...
    wavelet.internal.getmutexclopt(validopts,defaultopt,varargin);
useAdj = strcmpi(useAdj,"adjoint"); % Use a simple boolean value

% If original object size is not given, it is set to empty array.
if isempty(varargin)
    orgSize = [];
else
    orgSize = varargin{1};
end

if ~realA % Reorganize A into real valued array
    A = complex2cube(A);
end

% Default filters
params.Fsf = "nearsym5_7";
params.sf = 10;
params.sf = deblank(['qshift' num2str(params.sf)]);

% Get the level 1 filters
load(char(params.Fsf),'Lo*','Hi*');
% Get the q-shift filter
load(char(params.sf),'Lo*','Hi*');

if useAdj % Use adjoint in place of inverse
    % Biorthogonal reconstruction filters need to be replaced with
    % time-reversed decomposition filters. Since these filters are
    % symmetric, time reversal in not needed.
    LoR = LoD; 
    HiR = HiD;
end

% Switch to single if needed
if ~useDouble
    LoR = single(LoR);
    HiR = single(HiR);
    LoRa = single(LoRa);
    HiRa = single(HiRa);
end

% First level filters
g0 = LoR;
g1 = HiR;

% Levels >= 2
g0a = LoRa;
g1a = HiRa;

% Normalize analysis filters
gscale = 1 / norm(g0a,2);
% Tree A synthesis filters
g0a = g0a.* gscale;
g1a = g1a.* gscale;

% Tree B synthesis filters
g0b = g0a(end:-1:1);
g1b = g1a(end:-1:1);

% Debug cleanup
clear Hi* Lo*

while level > 1             % Go through levels in decreasing order
    if ~isempty(D{level-1}) % Get the desired size from lower level
        syh = size(D{level-1});
        prev_level_size = syh(1:4);
    else
        syh = size(D{level});
        prev_level_size = syh(1:4).*2;
    end
    % Obtain scaling coefficients through synthesis
    A = level2synthesis(A,D{level},g0a,g0b,g1a,g1b,prev_level_size);
    D{level} = [];   % Clear used detail coefficients to save memory.
    level = level-1;
end

% Special case where 1st level details were excluded
if level == 1 && isempty(D{1})
    rec = level1SynthNoHighpass(A,g0,orgSize);
end
% Synthesis from 1st level detail coefficients
if level == 1 && ~isempty(D{1})
    rec = level1SynthHighpass(A,D{1},g0,g1);
end
end

%------------------------------------------------------------------------
function Y = invOddEvenFilter(X,ha,hb,dim)
% Convolve one dimension of X using both filters. Interlace the values to
% double the length of the convolved direction.

% Case filters to column vectors
ha = ha(:);
hb = hb(:);
M = length(ha);
% The following will only work with even length filters
L = fix(M/2);

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

szX = size(X); % Input size

od = 2*szX(dim); % Operated dimension is doubled
szY = szX; szY(dim) = od;
Y = zeros(szY,class(X));

extIdx = wextend('ac','sym',(uint8(1:szX(dim))),L);
% Polyphase components of the filters
haOdd = ha(1:2:M);
haEven = ha(2:2:M);
hbOdd = hb(1:2:M);
hbEven = hb(2:2:M);

s = uint8(1:4:od); % Interlacing indicies

t = uint8(4:2:(szX(dim)+M)); % matIdx indicies
if dot(ha,hb) > 0
    ta = t; tb = t - 1;
else
    ta = t - 1; tb = t;
end
% Cell array to operate correct dimension
Ja = cell(1,4); Ja(:) = {':'};
Jb = Ja;
% Create cell arrays where correct indices are placed on dim'th place
Ia = Ja;
Ib = Ja;

if rem(L,2) == 1    % L is odd 
    ta = ta - 1; % Shited one to the left
    tb = tb - 1; % Shited one to the left
    
    % Switch "odd" and "even" polyphase filters around. Easiest way is to
    % use incorrect filter names!!!
    % Redo polyphase components of the filters
    haOdd = ha(2:2:M);
    haEven = ha(1:2:M);
    hbOdd = hb(2:2:M);
    hbEven = hb(1:2:M);
    
    % No need to shift ta and tb for "even" filters
    Ia{dim} = extIdx(ta);
    Ib{dim} = extIdx(tb);
else
    % ta and tb need to be shifted for "even" filters
    Ia{dim} = extIdx(ta-2);
    Ib{dim} = extIdx(tb-2);
end

% Even filter values
Ja{dim} = s+1;
Jb{dim} = s;
% Convolve
Y(Ja{:}) = convn(X(Ia{:}),haEven,'valid');
Y(Jb{:}) = convn(X(Ib{:}),hbEven,'valid');

% Odd filter values
Ja{dim} = s+3;
Jb{dim} = s+2;
Ia{dim} = extIdx(ta);
Ib{dim} = extIdx(tb);
% Convolve
Y(Ja{:}) = convn(X(Ia{:}),haOdd,'valid');
Y(Jb{:}) = convn(X(Ib{:}),hbOdd,'valid');
end

%------------------------------------------------------------------------
function Y = invOddEvenFilterOld(X,ha,hb,perm)
% Convolve one dimension of X using both filters. Interlace the values to
% double the length of the convolved direction.

% Permute X so that the first dimension is convolved
if ~isempty(perm); X = permute(X,perm); end

[r,c,s,t] = size(X); % Rows, Columns, Slices, Time steps
Y = zeros(2*r,c,s,t,class(X));

filtlen = length(ha);
% The following will only work with even length filters
L = fix(filtlen/2);
matIdx = wextend('ar','sym',(uint8(1:r))',L);
% Polyphase components of the filters
haOdd = ha(1:2:filtlen);
haEven = ha(2:2:filtlen);
hbOdd = hb(1:2:filtlen);
hbEven = hb(2:2:filtlen);

s = uint8(1:4:(r*2)); % Interlacing indicies

t = uint8(4:2:(r+filtlen)); % matIdx indicies
if dot(ha,hb) > 0
    ta = t; tb = t - 1;
else
    ta = t - 1; tb = t;
end
if rem(L,2) == 0    % L is even  
    
    % Convolve
    Y(s,:,:,:)   = convn(X(matIdx(tb-2),:,:,:),hbEven(:),'valid');
    Y(s+1,:,:,:) = convn(X(matIdx(ta-2),:,:,:),haEven(:),'valid');
    Y(s+2,:,:,:) = convn(X(matIdx(tb),:,:,:),hbOdd(:),'valid');
    Y(s+3,:,:,:) = convn(X(matIdx(ta),:,:,:),haOdd(:),'valid'); 
    
else % L is odd
    
    ta = ta - 1; % Shited one to the left
    tb = tb - 1; % Shited one to the left
    
    % Convolve (note the different order!)
    Y(s,:,:,:)   = convn(X(matIdx(tb),:,:,:),hbOdd(:),'valid');
    Y(s+1,:,:,:) = convn(X(matIdx(ta),:,:,:),haOdd(:),'valid');
    Y(s+2,:,:,:) = convn(X(matIdx(tb),:,:,:),hbEven(:),'valid');
    Y(s+3,:,:,:) = convn(X(matIdx(ta),:,:,:),haEven(:),'valid');
end
% Revert permutation
if ~isempty(perm); Y = ipermute(Y,perm); end
end

%-------------------------------------------------------------------------
function Z = complex2cube(Z)
% Reorganize the complex valued array into interlaced real valued array of
% twice the size.

global useAdj % Normalization differs for inverse and adjoint
if ~useAdj
    % Each orthant was divided by 2 in analysis, 8*1/2 = 4 hence 
    c = 0.25;
else
    % Adjoint must use same multiplier as analysis!
    c = 0.5;
end

sZ = size(Z);
% Split the array into real and imaginary parts of the 8 orthants
% Real parts
O1r = c*real(Z(:,:,:,:,1));
O2r = c*real(Z(:,:,:,:,2));
O3r = c*real(Z(:,:,:,:,3));
O4r = c*real(Z(:,:,:,:,4));
O5r = c*real(Z(:,:,:,:,5));
O6r = c*real(Z(:,:,:,:,6));
O7r = c*real(Z(:,:,:,:,7));
O8r = c*real(Z(:,:,:,:,8));
% Imaginary parts
O1i = c*imag(Z(:,:,:,:,1));
O2i = c*imag(Z(:,:,:,:,2));
O3i = c*imag(Z(:,:,:,:,3));
O4i = c*imag(Z(:,:,:,:,4));
O5i = c*imag(Z(:,:,:,:,5));
O6i = c*imag(Z(:,:,:,:,6));
O7i = c*imag(Z(:,:,:,:,7));
O8i = c*imag(Z(:,:,:,:,8));

% Allocate array for result
Z = zeros(2*sZ(1:4),class(O1r));

% Each orthant O_k is built from 8 functions Psi_l, 4 for real part and 4
% for imaginary. The sign of these functions is based on the tree-like
% structure and how the sign of the imaginary unit changes. To only obtain
% one function from the different orthants the signs have to match the
% specific function to obtain 8 copies of it and the other 7 functions
% cancel out.

% Combine real parts
Z(2:2:end,2:2:end,2:2:end,2:2:end) = +O1r+O2r+O3r+O4r+O5r+O6r+O7r+O8r; % Paaaa
Z(2:2:end,2:2:end,1:2:end,1:2:end) = -O1r-O2r-O3r-O4r+O5r+O6r+O7r+O8r; % Paabb
Z(2:2:end,1:2:end,2:2:end,1:2:end) = -O1r-O2r+O3r+O4r-O5r-O6r+O7r+O8r; % Pabab
Z(2:2:end,1:2:end,1:2:end,2:2:end) = -O1r-O2r+O3r+O4r+O5r+O6r-O7r-O8r; % Pabba
Z(1:2:end,2:2:end,2:2:end,1:2:end) = -O1r+O2r-O3r+O4r-O5r+O6r-O7r+O8r; % Pbaab
Z(1:2:end,2:2:end,1:2:end,2:2:end) = -O1r+O2r-O3r+O4r+O5r-O6r+O7r-O8r; % Pbaba
Z(1:2:end,1:2:end,2:2:end,2:2:end) = -O1r+O2r+O3r-O4r-O5r+O6r+O7r-O8r; % Pbbaa
Z(1:2:end,1:2:end,1:2:end,1:2:end) = +O1r-O2r-O3r+O4r-O5r+O6r+O7r-O8r; % Pbbbb

% Combine imaginary parts
Z(2:2:end,2:2:end,2:2:end,1:2:end) = +O1i+O2i+O3i+O4i+O5i+O6i+O7i+O8i; % Paaab
Z(2:2:end,2:2:end,1:2:end,2:2:end) = +O1i+O2i+O3i+O4i-O5i-O6i-O7i-O8i; % Paaba
Z(2:2:end,1:2:end,2:2:end,2:2:end) = +O1i+O2i-O3i-O4i+O5i+O6i-O7i-O8i; % Pabaa
Z(2:2:end,1:2:end,1:2:end,1:2:end) = -O1i-O2i+O3i+O4i+O5i+O6i-O7i-O8i; % Pabbb
Z(1:2:end,2:2:end,2:2:end,2:2:end) = +O1i-O2i+O3i-O4i+O5i-O6i+O7i-O8i; % Pbaaa
Z(1:2:end,2:2:end,1:2:end,1:2:end) = -O1i+O2i-O3i+O4i+O5i-O6i+O7i-O8i; % Pbabb
Z(1:2:end,1:2:end,2:2:end,1:2:end) = -O1i+O2i+O3i-O4i-O5i+O6i+O7i-O8i; % Pbbab
Z(1:2:end,1:2:end,1:2:end,2:2:end) = -O1i+O2i+O3i-O4i+O5i-O6i-O7i+O8i; % Pbbba
end

%-------------------------------------------------------------------------
function Y = level2synthesis(Yl,Yh,g0a,g0b,g1a,g1b,prev_level_size)
% From this we will only return the scaling coefficients at the next finer
% resolution level
LLLLsize = size(Yl);
% It's faster to work with two smaller arrays
Y1 = zeros([LLLLsize(1),2*LLLLsize(2:end)],class(Yh));
Y2 = Y1;
sY = size(Y1);

s2a = uint8(1:sY(2)/2);
s3a = uint8(1:sY(3)/2);
s4a = uint8(1:sY(4)/2);
s2b = s2a+uint8(sY(2)/2);
s3b = s3a+uint8(sY(3)/2);
s4b = s4a+uint8(sY(4)/2);

% Fill y array for synthesis
Y1(:,s2a,s3a,s4a) = Yl;

% There are 120 directions, 8 for each of the 15 filter setups
ind = reshape(uint8(1:120), [8,15]).'; % Helpful index matrix

Y2(:,s2a,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(1,:)));    % HLLL
Y1(:,s2b,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(2,:)));    % LHLL
Y2(:,s2b,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(3,:)));    % HHLL
Y1(:,s2a,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(4,:)));    % LLHL
Y2(:,s2a,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(5,:)));    % HLHL
Y1(:,s2b,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(6,:)));    % LHHL
Y2(:,s2b,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(7,:)));    % HHHL
Y1(:,s2a,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(8,:)));    % LLLH
Y2(:,s2a,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(9,:)));    % HLLH
Y1(:,s2b,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(10,:)));   % LHLH
Y2(:,s2b,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(11,:)));   % HHLH
Y1(:,s2a,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(12,:)));   % LLHH
Y2(:,s2a,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(13,:)));   % HLHH
Y1(:,s2b,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(14,:)));   % LHHH
Y2(:,s2b,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(15,:)));   % HHHH

size_curr_level = size(Yh);
size_curr_level = size_curr_level(1:4);

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T

% Filter dimension 1
perm = []; % Same as [1 2 3 4]
Yl = invOddEvenFilter(Y1,g0a,g0b,1);  % Lowpass
Yh = invOddEvenFilter(Y2,g1a,g1b,1);  % Highpass
Y = Yl + Yh;

% Filter dimension 2
perm = [2 1 3 4];
Yl = invOddEvenFilter(Y(:,s2a,:,:),g0a,g0b,2); % Lowpass
Yh = invOddEvenFilter(Y(:,s2b,:,:),g1a,g1b,2); % Highpass
Y = Yl + Yh;

% Filter dimension 3
perm = [3 2 1 4];
Yl = invOddEvenFilter(Y(:,:,s3a,:),g0a,g0b,3); % Lowpass
Yh = invOddEvenFilter(Y(:,:,s3b,:),g1a,g1b,3); % Highpass
Y = Yl + Yh;

% Filter dimension 4
perm = [4 2 3 1];
Yl = invOddEvenFilter(Y(:,:,:,s4a),g0a,g0b,4); % Lowpass
Yh = invOddEvenFilter(Y(:,:,:,s4b),g1a,g1b,4); % Highpass
Y = Yl + Yh;

% Now check if the size of the previous level is exactly twice the size
% of the current level. If it is exactly twice the size, the data was not
% extended at the previous level, if it is not, we have to remove the
% added row, column, and page dimensions.

% X
if prev_level_size(1) ~= 2*size_curr_level(1)
    Y = Y(2:end-1,:,:,:);
end
% Y
if  prev_level_size(2) ~= 2*size_curr_level(2)
    Y = Y(:,2:end-1,:,:);
end
% Z
if prev_level_size(3) ~= 2*size_curr_level(3)
    Y = Y(:,:,2:end-1,:);
end
% T
if prev_level_size(4) ~= 2*size_curr_level(4)
    Y = Y(:,:,:,2:end-1);
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
function Y = level1SynthNoHighpass(Y,g0,orgSize)
% Finest level synthesis using only scaling coefficients.

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T

% Filter dimension 1
perm = []; % Same as [1 2 3 4]
Y = oneDFilter(Y,g0,1);

% Filter dimension 2
perm = [2 1 3 4];
Y = oneDFilter(Y,g0,2);

% Filter dimension 3
perm = [3 2 1 4];
Y = oneDFilter(Y,g0,3);

% Filter dimension 4
perm = [4 2 3 1];
Y = oneDFilter(Y,g0,4);

% If the user specifies "excludeL1" at the input, it is possible that
% the output data size may not be correct. To correct for that, the user
% can provide the orignal data size as an input.

if ~isempty(orgSize)  % Original size is known  
    size_curr_level = size(Y);
    % Dimensions are cut if needed
    if orgSize(1) ~= size_curr_level(1)
        Y = Y(2:end-1,:,:,:);
    end
    if  orgSize(2) ~= size_curr_level(2)
        Y = Y(:,2:end-1,:,:);
    end
    if orgSize(3) ~= size_curr_level(3)
        Y = Y(:,:,2:end-1,:);
    end
    if orgSize(4) ~= size_curr_level(4)
        Y = Y(:,:,:,2:end-1);
    end
end
end

%-------------------------------------------------------------------------
function Y = level1SynthHighpass(Yl,Yh,g0,g1)
% Finest level synthesis using all coefficients.

LLLLsize = size(Yl);
% It's faster to work with two smaller arrays
Y1 = zeros([LLLLsize(1),2*LLLLsize(2:end)],class(Yh));
Y2 = Y1;

s2a = uint8(1:LLLLsize(2));
s3a = uint8(1:LLLLsize(3));
s4a = uint8(1:LLLLsize(4));
s2b = s2a+LLLLsize(2);
s3b = s3a+LLLLsize(3);
s4b = s4a+LLLLsize(4);

% Fill y array for synthesis
Y1(:,s2a,s3a,s4a) = Yl;
% There are 120 directions, 8 for each of the 15 filter setups
ind = reshape(uint8(1:120), [8,15]).'; % Helpful index matrix

Y2(:,s2a,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(1,:)));    % HLLL
Y1(:,s2b,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(2,:)));    % LHLL
Y2(:,s2b,s3a,s4a) = complex2cube(Yh(:,:,:,:,ind(3,:)));    % HHLL
Y1(:,s2a,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(4,:)));    % LLHL
Y2(:,s2a,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(5,:)));    % HLHL
Y1(:,s2b,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(6,:)));    % LHHL
Y2(:,s2b,s3b,s4a) = complex2cube(Yh(:,:,:,:,ind(7,:)));    % HHHL
Y1(:,s2a,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(8,:)));    % LLLH
Y2(:,s2a,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(9,:)));    % HLLH
Y1(:,s2b,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(10,:)));   % LHLH
Y2(:,s2b,s3a,s4b) = complex2cube(Yh(:,:,:,:,ind(11,:)));   % HHLH
Y1(:,s2a,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(12,:)));   % LLHH
Y2(:,s2a,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(13,:)));   % HLHH
Y1(:,s2b,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(14,:)));   % LHHH
Y2(:,s2b,s3b,s4b) = complex2cube(Yh(:,:,:,:,ind(15,:)));   % HHHH

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Notice that the size of Y# is halved at every step in the respective
% direction/dimension.

% Filter dimension 1
Yl = oneDFilter(Y1,g0,1);      % Lowpass
Yh = oneDFilter(Y2,g1,1);      % Highpass
Y = Yl + Yh;

% Filter dimension 2
Yl = oneDFilter(Y(:,s2a,:,:),g0,2);      % Lowpass
Yh = oneDFilter(Y(:,s2b,:,:),g1,2);      % Highpass
Y = Yl + Yh;

% Filter dimension 3
Yl = oneDFilter(Y(:,:,s3a,:),g0,3);      % Lowpass
Yh = oneDFilter(Y(:,:,s3b,:),g1,3);      % Highpass
Y = Yl + Yh;

% Filter dimension 4
Yl = oneDFilter(Y(:,:,:,s4a),g0,4);     % Lowpass
Yh = oneDFilter(Y(:,:,:,s4b),g1,4);     % Highpass
Y = Yl + Yh;
end
