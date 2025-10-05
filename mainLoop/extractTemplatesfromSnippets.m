function [wTEMP, wPCA] = extractTemplatesfromSnippets(rez, nPCs)
% this function is very similar to extractPCfromSnippets.
% outputs not just the PC waveforms, but also the template "prototype", 
% basically k-means clustering of 1D waveforms. 

ops = rez.ops;

if isfield(ops,'run_clean')
    run_clean = ops.run_clean;
else
    run_clean = false;
end

if isfield(ops,'max_learned_spikes')
    max_learned_spikes = ops.max_learned_spikes;
else
    max_learned_spikes = 1e5;
end

% skip every this many batches
nskip = getOr(ops, 'nskip', 25);

Nbatch      = rez.temp.Nbatch;
NT  	= ops.NT;
batchstart = 0:NT:NT*Nbatch;

fid = fopen(ops.fproc, 'r'); % open the preprocessed data file

k = 0;
dd = gpuArray.zeros(ops.nt0, 5e4, 'single'); % preallocate matrix to hold 1D spike snippets
for ibatch = 1:nskip:Nbatch
    offset = 2 * ops.Nchan*batchstart(ibatch);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [ops.Nchan NT], '*int16');
    dat = dat';

    % move data to GPU and scale it back to unit variance
    dataRAW = gpuArray(dat); % Is this necessary???
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;

    % 
    if run_clean
        electrode_map = ElectrodeMap(rez.xc, rez.yc,60.0,ops.whiteningRange);
        [row, col] = find_spike_templates(dat, electrode_map, ops);
    else
        % find isolated spikes from each batch
        [row, col] = isolated_peaks_new(dataRAW, ops);
    end


    

    % for each peak, get the voltage snippet from that channel
    clips = get_SpikeSample(dataRAW, row, col, ops, 0);

    c = sq(clips(:, :));

    if k+size(c,2)>size(dd,2)
        dd(:, 2*size(dd,2)) = 0;
    end

    dd(:, k + (1:size(c,2))) = c;
    k = k + size(c,2);
    if k > max_learned_spikes
        disp('Exceeded max iterations in extractTemplatesfromSnippets..')
        break;
    end
end
fclose(fid);

% discard empty samples
dd = dd(:, 1:k);

disp(['Template size: ',num2str(size(dd))])

% Make sure all of the spikes are negative at the peak.
if run_clean
    valid_spikes = logical(ones(1,size(dd,2)));
    min_window = 1;
    for ii = -min_window:min_window
        valid_spikes = valid_spikes & (dd(ops.nt0min+ii,:) < 0);
    end
    dd = dd(:, valid_spikes);
    % Normalize the spikes.
    dd = dd ./ sum(dd.^2,1).^.5;
end


dd = double(gather(dd));
[U, Sv, V] = svdecon(dd); % the PCs are just the left singular vectors of the waveforms

wPCA = gpuArray(single(U(:, 1:nPCs))); % take as many as needed
wPCA(:,1) = - wPCA(:,1) * sign(wPCA(ops.nt0min,1));  % adjust the arbitrary sign of the first PC so its negativity is downward

% Initialize the templates.
if run_clean
    dd = gpuArray(dd);
    % Project the spikes onto the PCs.
    projection = dd' * wPCA;
    % Sort the projections to get the initial templates.
    [~, idx] = sort(projection, 'descend');
    % Find the nPCs evenly spaced templates.
    idx_spaced = round(linspace(1, size(dd,2), nPCs));
    wTEMP = dd(:, idx(idx_spaced));
    wTEMP = wTEMP ./ sum(wTEMP.^2,1).^.5; % normalize them

    for i = 1:10
        % at each iteration, assign the waveform to its most correlated cluster
        cc = wTEMP' * dd;
        [amax, imax] = max(cc,[],1);
        for j = 1:nPCs
            wTEMP(:,j)  = dd(:,imax==j) * amax(imax==j)'; % weighted average to get new cluster means
        end
        wTEMP = wTEMP ./ sum(wTEMP.^2,1).^.5; % unit normalize
    end
    wTEMP = single(wTEMP);
else
    dd = gpuArray(single(dd));
    % initialize the template clustering with random waveforms
    % wTEMP = dd(:, randperm(size(dd,2), nPCs));
    wTEMP = dd(:, round(linspace(1, size(dd,2), nPCs)));
    wTEMP = wTEMP ./ sum(wTEMP.^2,1).^.5; % normalize them

    for i = 1:10
    % at each iteration, assign the waveform to its most correlated cluster
    cc = wTEMP' * dd;
    [amax, imax] = max(cc,[],1);
    for j = 1:nPCs
        wTEMP(:,j)  = dd(:,imax==j) * amax(imax==j)'; % weighted average to get new cluster means
    end
    wTEMP = wTEMP ./ sum(wTEMP.^2,1).^.5; % unit normalize
    end
end