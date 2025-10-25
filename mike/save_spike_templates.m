function results = save_spike_templates(template_dir, results)

% Define the field names for the ops and results structures.
ops_names = {'wTEMP','wPCA'};
results_names = {'simScore','iNeighPC','iNeigh','dWU','W','U','mu','Wraw','nsp'};

out_struct = struct();
for ii = 1:length(ops_names)
    tmp = results.ops.(ops_names{ii});
    if isgpuarray(tmp)
        tmp = gather(tmp); % Move data from GPU to CPU
    end
    out_struct.(ops_names{ii}) = tmp; % Initialize each field in out_struct
end

for ii = 1 : length(results_names)
    tmp = results.(results_names{ii});
    if isgpuarray(tmp)
        tmp = gather(tmp); % Move data from GPU to CPU
    end
    out_struct.(results_names{ii}) = tmp; % Initialize each field in out_struct
end

save(fullfile(template_dir,'templates.mat'),'-struct','out_struct')



