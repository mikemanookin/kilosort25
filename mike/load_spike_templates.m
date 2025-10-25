function results = load_spike_templates(template_dir, results)

% Define the field names for the ops and results structures.
ops_names = {'wTEMP','wPCA'};
results_names = {'simScore','iNeighPC','iNeigh','dWU','W','U','mu','Wraw','nsp'};

% Load the templates.
template_filepath = fullfile(template_dir,'templates.mat');
if exist(template_filepath, "file")
    out_struct = load(template_filepath);
    for ii = 1:length(ops_names)
         results.ops.(ops_names{ii}) = out_struct.(ops_names{ii}); 
    end
    for ii = 1 : length(results_names)
         results.(results_names{ii}) = out_struct.(results_names{ii}); % Initialize each field in out_struct
    end
else
    error(["File does not exist at path: ", template_filepath]);
end

