function [row, col, mu] = isolated_peaks_new(S1, ops)
% takes a matrix of timepoints by channels S1
% outputs threshold crossings that are relatively isolated from other peaks
% outputs row, column and magnitude of the threshold crossing
loc_range = getOr(ops, 'loc_range', [5 4]);
long_range = getOr(ops, 'long_range', [30 6]);
Th = ops.spkTh;
nt0 = ops.nt0;

% loc_range = [3  1];
% long_range = [30  6];

% finding the local minimum in a sliding window within plus/minus loc_range extent
% across time and across channels
if sum(loc_range) > 2
    smin = my_min(S1, loc_range, [1 2]); % This should be a logical matrix.
    peaks = (S1<smin+1e-3 & S1<Th); % the peaks are samples that achieve this local minimum, AND have negativities less than a preset threshold
else % Local minimum calculation.
    local_min = find_local_minima(S1);
    peaks = local_min & (S1 < Th);
end

% Logical matrix of detected spikes
% peaks = (S1<smin+1e-3 & S1<Th);

% only take local peaks that are isolated from other local peaks
if sum(long_range) > 2 && sum(loc_range) > 2
    sum_peaks = my_sum(peaks, long_range, [1 2]); % if there is another local peak close by, this sum will be at least 2
    % peaks = single(peaks) .* (sum_peaks<1.2) .* S1; % set to 0 peaks that are not isolated, and multiply with the voltage values
    peaks = peaks & (sum_peaks<1.2);
end

% exclude temporal buffers
peaks([1:nt0 end-nt0:end], :) = 0;

[row, col] = find(peaks == 1); % find the non-zero peaks, and take their amplitudes

if nargout > 2
    % Get the amplitudes at the peak indices.
    mu = S1(peaks == 1);
    mu = -mu; % invert the sign of the amplitudes
end
