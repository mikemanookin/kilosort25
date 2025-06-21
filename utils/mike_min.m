function local_minima = mike_min(this_data)
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

dt = zeros(size(this_data));
for ii=1:2
    idx = 1:size(dt,1)-ii;
    dt(idx,:) = this_data(ii+1:end,:) - this_data(1:end-ii,:);
end
dt = dt/loc_range(1);

% Take the second derivative.
local_minima = zeros(size(this_data));
% ddt(1:end-1,:) = diff(dt > 0,1,1);
local_minima(2:end,:) = (dt(1:end-1,:) <= 0) & (dt(2:end,:) >= 0);


