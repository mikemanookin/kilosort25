function S1 = my_sum(S1, sig, varargin)
%MY_SUM Apply a running sum filter across specified dimensions.
%
%   S1_OUT = MY_SUM(S1, SIG) applies a 1D running sum filter to the input
%   array S1 along the default dimension (dim = 2), using a filter width
%   of 2*SIG + 1. SIG specifies the half-width of the summation window.
%
%   S1_OUT = MY_SUM(S1, SIG, DIMS) applies the filter along specified
%   dimensions. DIMS is a vector indicating the axes to be filtered.
%
%   SIG can be:
%     - A scalar: same filter width is applied along all specified axes
%     - A vector: each element specifies the filter width for the corresponding axis
%
%   Example:
%     % Apply a running sum with half-width 2 along columns
%     A = rand(10, 100);
%     B = my_sum(A, 2);
%
%     % Apply different filter widths along dimensions 1 and 2
%     C = rand(20, 50, 30);
%     D = my_sum(C, [3, 2], [1, 2]);
%
%   Notes:
%   - The function pads the array with zeros to avoid edge artifacts.
%   - Useful for computing local totals or pre-smoothing data before other
%     operations like thresholding or normalization.
%
%   Inputs:
%     S1    - Input matrix or n-dimensional array
%     sig   - Scalar or vector of half-window sizes (number of bins on either side)
%     varargin{1} (optional) - Vector of dimensions to apply filtering
%
%   Output:
%     S1    - Output array with running sum applied

% Set default filtering dimension to 2
idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end

% Determine filter widths per dimension
if numel(idims) > 1 && numel(sig) > 1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

% Apply running sum along each specified dimension
for i = 1:length(idims)
    sig = sigall(i);          % Half-width for current dimension
    idim = idims(i);          % Current dimension index
    Nd = ndims(S1);           % Number of dimensions in input

    % Bring the dimension to the front for filtering
    S1 = permute(S1, [idim, 1:idim-1, idim+1:Nd]);

    % Store original reshaped size
    dsnew = size(S1);

    % Collapse all other dimensions for 2D operation
    S1 = reshape(S1, size(S1, 1), []);
    dsnew2 = size(S1);

    % Pad with zeros on top and bottom
    S1 = cat(1, zeros([sig, dsnew2(2)]), S1, zeros([sig, dsnew2(2)]));

    % Initialize the sum accumulator
    Ssum = S1(1:dsnew2(1), :);
    for j = 1:2*sig
        Ssum = Ssum + S1(j + (1:dsnew2(1)), :);
    end

    % Reshape to original dimensionality
    S1 = reshape(Ssum, dsnew);

    % Restore original dimension order
    S1 = permute(S1, [2:idim, 1, idim+1:Nd]);
end
