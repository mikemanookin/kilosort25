classdef ElectrodeMap < handle
    properties
        x_coords % Property to store x coordinates of electrodes
        y_coords % Property to store y coordinates of electrodes
        electrode_spacing % Property to store spacing between electrodes
        num_nearest_channels % Property to store number of nearest channels
        nearest_channels % Property to store array of nearest channels.
        one_away_channels % Property to store channels located one electrode away.
        number_of_channels % Property to store number of electrodes/channels
    end

    methods
        function obj = ElectrodeMap(x_coords, y_coords, electrode_spacing, num_nearest_channels)
            obj.x_coords = x_coords;
            obj.y_coords = y_coords;
            obj.electrode_spacing = electrode_spacing;
            obj.num_nearest_channels = num_nearest_channels;
            obj.number_of_channels = length(x_coords);
            obj.get_one_away_channels();
        end

        function get_one_away_channels(obj)
            obj.one_away_channels = nan(obj.number_of_channels, 7);
            for ii = 1 : obj.number_of_channels
                % Calculate the distance between the channels and target.
                distances = sqrt((obj.x_coords(ii) - obj.x_coords).^2 + (obj.y_coords(ii) - obj.y_coords).^2);
                [sorted_distances, sortedIndices] = sort(distances);
                sortedIndices = sortedIndices(sorted_distances < obj.electrode_spacing*1.25);
                obj.one_away_channels(ii, 1:length(sortedIndices)) = sortedIndices'; % Store indices of channels one electrode away
            end
        end

        function get_nearest_channels(obj)
        end
    end
end


