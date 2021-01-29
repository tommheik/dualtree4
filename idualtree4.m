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
%   This code is heavily based on the IDUALTREE3-function.
%
%   Tommi Heikkilä
%   Created 15.5.2020
%   Last edited 25.1.2021

% A should be a 4-D matrix output from dualtree3
validateattributes(A,{'numeric'},{'real','nonempty','finite','ndims',4},...
    'IDUALTREE4','A');
% D should be a cell array of length at least one (two)
validateattributes(D,{'cell'},{'nonempty'},'IDUALTREE4','D');

% Obtain the level of the transform
level = length(D);

% Use single precision if needed
useDouble = 1;
if isa(A, 'single')
    useDouble = 0;
end

% If original object size is not given, it is set to empty array.
if nargin < 3
    orgSize = [];
else
    orgSize = varargin{1};
end

% Default filters
params.Fsf = "nearsym5_7";
params.sf = 10;
params.sf = deblank(['qshift' num2str(params.sf)]);

% Get the level 1 filters
load(char(params.Fsf));
% Get the q-shift filter
load(char(params.sf));

% Switch to single if needed
if ~useDouble
    LoR = single(LoR);
    HiR = single(HiR);
    LoRa = single(LoRa);
    HiRa = single(HiRa);
end

% First level filters
g0o = LoR;
g1o = HiR;

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

if level>=2
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
end

% Special case where 1st level details were excluded
if level == 1 && isempty(D{1})
    A = level1SynthNoHighpass(A,g0o,orgSize);
end
% Synthesis from 1st level detail coefficients
if level == 1 && ~isempty(D{1})
    A = level1SynthHighpass(A,D{1},g0o,g1o);
end

% Return the recomposed 4-D object
rec = A;
end

%------------------------------------------------------------------------
function yrec = invColumnFilter(x,ha,hb)
[r,c] = size(x);
yrec = zeros(2*r,c);
filtlen = length(ha);
% The following will just work with even length filters
L = fix(filtlen/2);
matIdx = wextend('ar','sym',(uint8(1:r))',L);
% Polyphase components of the filters
haOdd = ha(1:2:filtlen);
haEven = ha(2:2:filtlen);
hbOdd = hb(1:2:filtlen);
hbEven = hb(2:2:filtlen);
s = uint8(1:4:(r*2));
if rem(L,2) == 0
    
    t = uint8(4:2:(r+filtlen));
    if dot(ha,hb) > 0
        ta = t; tb = t - 1;
    else
        ta = t - 1; tb = t;
    end
    
    yrec(s,:)   = conv2(x(matIdx(tb-2),:),haEven(:),'valid');
    yrec(s+1,:) = conv2(x(matIdx(ta-2),:),hbEven(:),'valid');
    yrec(s+2,:) = conv2(x(matIdx(tb),:),haOdd(:),'valid');
    yrec(s+3,:) = conv2(x(matIdx(ta),:),hbOdd(:),'valid'); 
    
else
    
    t = uint8(3:2:(r+filtlen-1));
    if dot(ha,hb) > 0
        ta = t; tb = t - 1;
    else
        ta = t - 1; tb = t;
    end
        
    s = uint8(1:4:(r*2));

    yrec(s,:)   = conv2(x(matIdx(tb),:),haOdd(:),'valid');
    yrec(s+1,:) = conv2(x(matIdx(ta),:),hbOdd(:),'valid');
    yrec(s+2,:) = conv2(x(matIdx(tb),:),haEven(:),'valid');
    yrec(s+3,:) = conv2(x(matIdx(ta),:),hbEven(:),'valid');
end
end

%-------------------------------------------------------------------------
function y = complex2cube(z)

sz = size(z);
% Split the array into real and imaginary parts of the 8 orthants
% Real parts
O1r = real(z(:,:,:,:,1));
O2r = real(z(:,:,:,:,2));
O3r = real(z(:,:,:,:,3));
O4r = real(z(:,:,:,:,4));
O5r = real(z(:,:,:,:,5));
O6r = real(z(:,:,:,:,6));
O7r = real(z(:,:,:,:,7));
O8r = real(z(:,:,:,:,8));
% Imaginary parts
O1i = imag(z(:,:,:,:,1));
O2i = imag(z(:,:,:,:,2));
O3i = imag(z(:,:,:,:,3));
O4i = imag(z(:,:,:,:,4));
O5i = imag(z(:,:,:,:,5));
O6i = imag(z(:,:,:,:,6));
O7i = imag(z(:,:,:,:,7));
O8i = imag(z(:,:,:,:,8));

% Save memory
clear z

% Allocate array for result
y = zeros(2*sz(1:4));

% Each orthant O_k is built from 8 functions Psi_l, 4 for real part and 4
% for imaginary. The sign of these functions is based on the tree-like
% structure and how the sign of the imaginary unit changes. To only obtain
% one function from the different orthants the signs have to match the
% specific function to obtain 8 copies of it and the other 7 functions
% cancel out.

% Combine real parts
y(2:2:end,2:2:end,2:2:end,2:2:end) = +O1r+O2r+O3r+O4r+O5r+O6r+O7r+O8r;
y(2:2:end,2:2:end,1:2:end,1:2:end) = -O1r-O2r-O3r-O4r+O5r+O6r+O7r+O8r;
y(2:2:end,1:2:end,2:2:end,1:2:end) = -O1r-O2r+O3r+O4r-O5r-O6r+O7r+O8r;
y(2:2:end,1:2:end,1:2:end,2:2:end) = -O1r-O2r+O3r+O4r+O5r+O6r-O7r-O8r;
y(1:2:end,2:2:end,2:2:end,1:2:end) = -O1r+O2r-O3r+O4r-O5r+O6r-O7r+O8r;
y(1:2:end,2:2:end,1:2:end,2:2:end) = -O1r+O2r-O3r+O4r+O5r-O6r+O7r-O8r;
y(1:2:end,1:2:end,2:2:end,2:2:end) = -O1r+O2r+O3r-O4r-O5r+O6r+O7r-O8r;
y(1:2:end,1:2:end,1:2:end,1:2:end) = +O1r-O2r-O3r+O4r-O5r+O6r+O7r-O8r;

% Combine imaginary parts
y(2:2:end,2:2:end,2:2:end,1:2:end) = +O1i+O2i+O3i+O4i+O5i+O6i+O7i+O8i;
y(2:2:end,2:2:end,1:2:end,2:2:end) = +O1i+O2i+O3i+O4i-O5i-O6i-O7i-O8i;
y(2:2:end,1:2:end,2:2:end,2:2:end) = +O1i+O2i-O3i-O4i+O5i+O6i-O7i-O8i;
y(2:2:end,1:2:end,1:2:end,1:2:end) = -O1i-O2i+O3i+O4i+O5i+O6i-O7i-O8i;
y(1:2:end,2:2:end,2:2:end,2:2:end) = +O1i-O2i+O3i-O4i+O5i-O6i+O7i-O8i;
y(1:2:end,2:2:end,1:2:end,1:2:end) = -O1i+O2i-O3i+O4i+O5i-O6i+O7i-O8i;
y(1:2:end,1:2:end,2:2:end,1:2:end) = -O1i+O2i+O3i-O4i-O5i+O6i+O7i-O8i;
y(1:2:end,1:2:end,1:2:end,2:2:end) = -O1i+O2i+O3i-O4i+O5i-O6i-O7i+O8i;

% Each orthant was divided by 2 in analysis, 8*1/2 = 4 hence 
y = 0.25*y;
end

%-------------------------------------------------------------------------
function y = level2synthesis(yl,yh,g0a,g0b,g1a,g1b,prev_level_size)
% From this we will only return the scaling coefficients at the next finer
% resolution level
LLLLsize = size(yl);
if isa(yl,'single') || isa(yh,'single')
    yl = single(yl);
    yh = single(yh);
    y = zeros(2*LLLLsize,'single');
else
    y = zeros(2*LLLLsize);
end
sy = size(y);

s1a = uint8(1:sy(1)/2);
s2a = uint8(1:sy(2)/2);
s3a = uint8(1:sy(3)/2);
s4a = uint8(1:sy(4)/2);
s1b = s1a+uint8(sy(1)/2);
s2b = s2a+uint8(sy(2)/2);
s3b = s3a+uint8(sy(3)/2);
s4b = s4a+uint8(sy(4)/2);

% Fill y array for synthesis
y(s1a,s2a,s3a,s4a) = yl;
clear yl
% There are 120 directions, 8 for each of the 15 filter setups
ind = reshape(uint8(1:120), [8,15]).'; % Helpful index matrix

y(s1b,s2a,s3a,s4a) = complex2cube(yh(:,:,:,:,ind(1,:)));    % HLLL
y(s1a,s2b,s3a,s4a) = complex2cube(yh(:,:,:,:,ind(2,:)));    % LHLL
y(s1b,s2b,s3a,s4a) = complex2cube(yh(:,:,:,:,ind(3,:)));    % HHLL
y(s1a,s2a,s3b,s4a) = complex2cube(yh(:,:,:,:,ind(4,:)));    % LLHL
y(s1b,s2a,s3b,s4a) = complex2cube(yh(:,:,:,:,ind(5,:)));    % HLHL
y(s1a,s2b,s3b,s4a) = complex2cube(yh(:,:,:,:,ind(6,:)));    % LHHL
y(s1b,s2b,s3b,s4a) = complex2cube(yh(:,:,:,:,ind(7,:)));    % HHHL
y(s1a,s2a,s3a,s4b) = complex2cube(yh(:,:,:,:,ind(8,:)));    % LLLH
y(s1b,s2a,s3a,s4b) = complex2cube(yh(:,:,:,:,ind(9,:)));    % HLLH
y(s1a,s2b,s3a,s4b) = complex2cube(yh(:,:,:,:,ind(10,:)));   % LHLH
y(s1b,s2b,s3a,s4b) = complex2cube(yh(:,:,:,:,ind(11,:)));   % HHLH
y(s1a,s2a,s3b,s4b) = complex2cube(yh(:,:,:,:,ind(12,:)));   % LLHH
y(s1b,s2a,s3b,s4b) = complex2cube(yh(:,:,:,:,ind(13,:)));   % HLHH
y(s1a,s2b,s3b,s4b) = complex2cube(yh(:,:,:,:,ind(14,:)));   % LHHH
y(s1b,s2b,s3b,s4b) = complex2cube(yh(:,:,:,:,ind(15,:)));   % HHHH

% yh no longer needed so clear it to save memory!
size_curr_level = size(yh);
size_curr_level = size_curr_level(1:4);
clear yh

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Filter dimensions 1 and 2
% Note that we loop through full number of slices and time steps!
for timeidx = 1:sy(4)
    for sliceidx = 1:sy(3)
        % 1st dimension
        ytmp = invColumnFilter(y(s1a,:,sliceidx,timeidx),g0b,g0a) + ...
            invColumnFilter(y(s1b,:,sliceidx,timeidx),g1b,g1a);
        % 2nd dimension
        % ytmp has to be transposed before and after filtering
        y(:,:,sliceidx,timeidx) = invColumnFilter(ytmp(:,s2a).',g0b,g0a).' + ...
            invColumnFilter(ytmp(:,s2b).',g1b,g1a).';
    end
end

% Filter dimensions 3 and 4
% Note that we loop through  full number of rows and columns!
for colidx = 1:sy(2)
    for rowidx = 1:sy(1)
        % 3rd dimension
        ytmp = invColumnFilter(squeeze(y(rowidx,colidx,s3a,:)),g0b,g0a) + ...
            invColumnFilter(squeeze(y(rowidx,colidx,s3b,:)),g1b,g1a);
        % 4th dimension
        % ytmp has to be transposed before and after filtering
        y(rowidx,colidx,:,:) = reshape((invColumnFilter(ytmp(:,s4a).',g0b,g0a) + ...
            invColumnFilter(ytmp(:,s4b).',g1b,g1a)).',[1,1,sy(3),sy(4)]);
    end
end
clear ytmp

% Now check if the size of the previous level is exactly twice the size
% of the current level. If it is exactly twice the size, the data was not
% extended at the previous level, if it is not, we have to remove the
% added row, column, and page dimensions.

% X
if prev_level_size(1) ~= 2*size_curr_level(1)
    y = y(2:end-1,:,:,:);
end
% Y
if  prev_level_size(2) ~= 2*size_curr_level(2)
    y = y(:,2:end-1,:,:);
end
% Z
if prev_level_size(3) ~= 2*size_curr_level(3)
    y = y(:,:,2:end-1,:);
end
% T
if prev_level_size(4) ~= 2*size_curr_level(4)
    y = y(:,:,:,2:end-1);
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
function yl = level1SynthNoHighpass(yl,g0o,orgSize)
% We use no highpass filter here
LLLLsize = size(yl);
s1 = uint8(1:LLLLsize(1));
s2 = uint8(1:LLLLsize(2));
s3 = uint8(1:LLLLsize(3));
s4 = uint8(1:LLLLsize(4));

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Filter dimensions 1 and 2
for sliceidx = s3
    for timeidx = s4
        % 1st dimension
        y = columnFilter(yl(:,:,sliceidx,timeidx),g0o);
        % 2nd dimension
        % y is tranposed before and after filtering.
        yl(:,:,sliceidx,timeidx) = columnFilter(y.',g0o).';
    end
end

% Filter dimensions 3 and 4
for rowidx = s1
    for colidx = s2
        % 3rd dimension
        y = columnFilter(squeeze(yl(rowidx,colidx,:,:)),g0o);
        % 4th dimension
        % y is transposed before and after filtering and then extended to
        % a part of 4-D array.
        yl(rowidx,colidx,:,:) = reshape(columnFilter(y.',g0o).',[1,1,LLLLsize(3),LLLLsize(4)]);
    end
end

% If the user specifies "excludeL1" at the input, it is possible that
% the output data size may not be correct. To correct for that, the user
% can provide the orignal data size as an input.

if ~isempty(orgSize)  % Original size is known  
    size_curr_level = size(yl);
    % Dimensions are cut if needed
    if orgSize(1) ~= size_curr_level(1)
        yl = yl(2:end-1,:,:,:);
    end
    if  orgSize(2) ~= size_curr_level(2)
        yl = yl(:,2:end-1,:,:);
    end
    if orgSize(3) ~= size_curr_level(3)
        yl = yl(:,:,2:end-1,:);
    end
    if orgSize(4) ~= size_curr_level(4)
        yl = yl(:,:,:,2:end-1);
    end
end
end

%-------------------------------------------------------------------------
function y = level1SynthHighpass(yl,yh,g0o,g1o)

LLLLsize = size(yl);
if isa(yl,'single') || isa(yh,'single')
    yl = single(yl);
    yh = single(yh);
    y = zeros(2*LLLLsize,'single');
else
    y = zeros(2*LLLLsize);
end
sy = size(y);

x1a = uint8(1:LLLLsize(1));
x2a = uint8(1:LLLLsize(2));
x3a = uint8(1:LLLLsize(3));
x4a = uint8(1:LLLLsize(4));
x1b = x1a+LLLLsize(1);
x2b = x2a+LLLLsize(2);
x3b = x3a+LLLLsize(3);
x4b = x4a+LLLLsize(4);

% Fill y array for synthesis
y(x1a,x2a,x3a,x4a) = yl;
clear yl
% There are 120 directions, 8 for each of the 15 filter setups
ind = reshape(uint8(1:120), [8,15]).'; % Helpful index matrix

y(x1b,x2a,x3a,x4a) = complex2cube(yh(:,:,:,:,ind(1,:)));    % HLLL
y(x1a,x2b,x3a,x4a) = complex2cube(yh(:,:,:,:,ind(2,:)));    % LHLL
y(x1b,x2b,x3a,x4a) = complex2cube(yh(:,:,:,:,ind(3,:)));    % HHLL
y(x1a,x2a,x3b,x4a) = complex2cube(yh(:,:,:,:,ind(4,:)));    % LLHL
y(x1b,x2a,x3b,x4a) = complex2cube(yh(:,:,:,:,ind(5,:)));    % HLHL
y(x1a,x2b,x3b,x4a) = complex2cube(yh(:,:,:,:,ind(6,:)));    % LHHL
y(x1b,x2b,x3b,x4a) = complex2cube(yh(:,:,:,:,ind(7,:)));    % HHHL
y(x1a,x2a,x3a,x4b) = complex2cube(yh(:,:,:,:,ind(8,:)));    % LLLH
y(x1b,x2a,x3a,x4b) = complex2cube(yh(:,:,:,:,ind(9,:)));    % HLLH
y(x1a,x2b,x3a,x4b) = complex2cube(yh(:,:,:,:,ind(10,:)));   % LHLH
y(x1b,x2b,x3a,x4b) = complex2cube(yh(:,:,:,:,ind(11,:)));   % HHLH
y(x1a,x2a,x3b,x4b) = complex2cube(yh(:,:,:,:,ind(12,:)));   % LLHH
y(x1b,x2a,x3b,x4b) = complex2cube(yh(:,:,:,:,ind(13,:)));   % HLHH
y(x1a,x2b,x3b,x4b) = complex2cube(yh(:,:,:,:,ind(14,:)));   % LHHH
y(x1b,x2b,x3b,x4b) = complex2cube(yh(:,:,:,:,ind(15,:)));   % HHHH
clear yh

% Filter dimensions in reverse order compared to analysis: X->Y->Z->T
% Filter dimensions 1 and 2
% Note that we loop through doubled number of slices and time steps!
for timeidx = 1:sy(4)
    for sliceidx = 1:sy(3)
        % 1st dimension
        ytmp = columnFilter(y(x1a,:,sliceidx,timeidx),g0o) + ...
            columnFilter(y(x1b,:,sliceidx,timeidx),g1o);
        % 2nd dimension
        % ytmp has to be transposed before and after filtering
        y(x1a,x2a,sliceidx,timeidx) = columnFilter(ytmp(:,x2a).',g0o).' + ...
            columnFilter(ytmp(:,x2b).',g1o).';
    end
end
        
% Filter dimensions 3 and 4
% Note that we loop through  normal number of rows and columns!
for colidx = 1:LLLLsize(2)
    for rowidx = 1:LLLLsize(1)
        % 3rd dimension
        ytmp = columnFilter(squeeze(y(rowidx,colidx,x3a,:)),g0o) + ...
            columnFilter(squeeze(y(rowidx,colidx,x3b,:)),g1o);
        % 4th dimension
        % ytmp has to be transposed before and after filtering
        y(rowidx,colidx,x3a,x4a) = reshape((columnFilter(ytmp(:,x4a).',g0o) + ...
            columnFilter(ytmp(:,x4b).',g1o)).',[1,1,LLLLsize(3),LLLLsize(4)]);
    end
end

% Return synthezised object
y = y(x1a,x2a,x3a,x4a);
end
