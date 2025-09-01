function dataRAW = clip_out_spikes(dataRAW, spike_row, spike_col, spike_amp, varargin)

ip = inputParser();
ip.addParameter('sample_rate', 20000, @(x)isnumeric(x)); % Spike threshold
ip.parse(varargin{:});

sample_rate = ip.Results.sample_rate;


% Sort spikes by amplitude in decending order.
[~, spike_idx] = sort(spike_amp,'ascend');
spike_row = spike_row(spike_idx);
spike_col = spike_col(spike_idx);
spike_amp = spike_amp(spike_idx);

%% Clip out the spikes.
half_width = 21;
% data_copy = dataRAW;
% [num_rows, ~] = size(dataRAW);
% 
% for ii = 1:length(spike_row)
%     r = spike_row(ii);
%     c = spike_col(ii);
% 
%     % Define window bounds
%     c_start = max(r - half_width, 1);
%     c_end = min(r + half_width, num_rows);
% 
%     % Compute average of edge values
%     % window_values = dataRAW([c_start, c_end], c);
%     window_values = data_copy([c_start, c_end], c);
%     % replacement_value = mean(window_values);
% 
%     % Linear interpolation.
%     % window_data = dataRAW(c_start:c_end, c);
%     window_slope = diff(window_values)/(c_end-c_start);
%     replacement_value = (0:(c_end-c_start))*window_slope + window_values(1);
% 
%     % Replace values in the window
%     data_copy(c_start:c_end, c) = replacement_value;
% end

data_copy = clip_spikes(dataRAW, spike_row, spike_col, half_width);

[b1, a1] = butter(3, 500/sample_rate*2, 'low');
data_low = filter(b1, a1, data_copy); % causal forward filter
% datr = filter(b1, a1, d_filt.data); % causal forward filter
data_low = flipud(data_low); % reverse time
data_low = filter(b1, a1, data_low); % causal forward filter again
data_low = flipud(data_low); % reverse time back

% Subtract the low-frequency drift.
dataRAW = dataRAW - data_low;

% figure(1); clf;
% hold on
% plot(dataRAW(:,target_channel));
% plot(data_copy(:,target_channel));
% plot(data_low(:,target_channel),'k','LineWidth',1);
% hold off;
% set(gca,'XLim',[5800,6200])