classdef TemplateFinder < handle
    properties
        Templates % Property to store template data
    end

    methods
        function obj = TemplateFinder(templates)
            obj.Templates = templates;
        end

        function get_nearest_channels(x_coords, y_corrds, search_radius)
        end
    end
end


