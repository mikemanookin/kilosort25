function local_minima = mike_min(this_data)


% Calculate the first derivative.
dt = zeros(size(this_data));
dt(1:size(dt,1)-1,:) = this_data(2:end,:) - this_data(1:end-1,:);
% // for ii=1:2
% //     idx = 1:size(dt,1)-ii;
% //     dt(idx,:) = this_data(ii+1:end,:) - this_data(1:end-ii,:);
% // end
% // dt = dt/2;

% Take the second derivative.
local_minima = zeros(size(this_data));
local_minima(2:end,:) = (dt(1:end-1,:) <= 0) & (dt(2:end,:) >= 0);




