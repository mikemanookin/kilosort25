function [spike_row, spike_col, spike_amp] = find_spikes(dataRAW, varargin)

ip = inputParser();
ip.addParameter('threshold', -100, @(x)isfloat(x)); % Spike threshold
ip.parse(varargin{:});

threshold = ip.Results.threshold;

% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

% loc_range = [3,1];
% smin = my_min(dataRAW, loc_range, [1 2]);

spikes = find_local_minima(dataRAW) & (dataRAW < threshold);

% Find the row, column where the spikes live.
[spike_row, spike_col] = find(spikes);

if nargout > 2
    spike_amp = dataRAW(spikes);
end