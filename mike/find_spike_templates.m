function [row, col] = find_spike_templates(dat, electrode_map, ops)

Th = ops.spkTh;
nt0 = ops.nt0;

local_min = find_local_minima(dat);
peaks = local_min & (dat < Th);

% Loop through the peaks and take only values that are local minima.
for ii = 1 : size(dat,2)
    if sum(peaks(:,ii)) > 0
        neighbors = electrode_map.one_away_channels(ii,:);
        neighbors = neighbors(~isnan(neighbors));
        neighbor_min = min(dat(:,neighbors),[],2);
        peaks(:, ii) = peaks(:, ii) & (dat(:, ii) <= neighbor_min);
    end
end

% exclude temporal buffers
peaks([1:nt0 end-nt0:end], :) = 0;

[row, col] = find(peaks == 1); % find the non-zero peaks, and take their amplitudes
