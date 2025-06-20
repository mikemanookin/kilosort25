function S1 = my_min(S1, sig, varargin)
%MY_MIN Apply a running minimum filter across specified dimensions.
%
%   S1_OUT = MY_MIN(S1, SIG) applies a 1D running minimum filter to the
%   input matrix S1 along the default dimension (dim = 2) using a
%   window size of 2*SIG + 1. SIG specifies the half-width of the filter
%   window.
%
%   S1_OUT = MY_MIN(S1, SIG, DIMS) allows filtering along specific
%   dimensions. DIMS is a vector specifying which axes to filter.
%
%   SIG can be a scalar or a vector:
%     - If SIG is a scalar, the same filter size is applied across all
%       specified dimensions.
%     - If SIG is a vector with the same length as DIMS, each dimension
%       will be filtered with its corresponding SIG value.
%
%   Example:
%     % Apply a running minimum with half-width 2 along columns
%     A = rand(10, 100);
%     B = my_min(A, 2);
%
%     % Apply different filter sizes along dimensions 1 and 2
%     C = rand(20, 50, 30);
%     D = my_min(C, [3, 2], [1, 2]);
%
%   Notes:
%   - Padding with Inf ensures that values outside the filter window do not
%     influence the result.
%   - The function uses reshaping and permutation to apply the filter
%     along arbitrary dimensions.

% Default to filtering along dimension 2 if no dimension is specified
idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end

% Handle scalar or vector SIG specification
if numel(idims)>1 && numel(sig)>1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

% Loop through each specified dimension
for i = 1:length(idims)
    sig = sigall(i); % Filter half-width for current dimension
    idim = idims(i); % Current dimension to filter
    Nd = ndims(S1); % Number of dimensions in input

    % Move target dimension to the front
    S1 = permute(S1, [idim, 1:idim-1, idim+1:Nd]);

    dsnew = size(S1);

    % Reshape to 2D for filtering: (target dim) x (everything else)
    S1 = reshape(S1, size(S1,1), []);
    dsnew2 = size(S1);

    % Pad with Inf to avoid edge effects
    S1 = cat(1, Inf*ones([sig, dsnew2(2)]),S1, Inf*ones([sig, dsnew2(2)]));
    
    % Apply running minimum filter
    Smax = S1(1:dsnew2(1), :);
    for j = 1:2*sig
        Smax = min(Smax, S1(j + (1:dsnew2(1)), :));
    end

    % Restore original shape
    S1 = reshape(Smax, dsnew);

    % Permute dimensions back to original order
    S1 = permute(S1, [2:idim, 1, idim+1:Nd]);
end
