function pre_filter_data()

% ip = inputParser();
% ip.addParameter('threshold', -100, @(x)isfloat(x)); % Spike threshold
% ip.parse(varargin{:});
% 
% threshold = ip.Results.threshold;

in_file_path = '/data/data/sorted/20250306C/doves_images_bak.bin';
out_file_path = '/data/data/sorted/20250306C/doves_images.bin';

NT       = 65472; %ops.NT ; % number of timepoints per batch
NchanTOT = 512; %ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc
NTbuff = 65664;
bytes       = get_file_size(in_file_path);
nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints

Nbatch      = ceil(nTimepoints /NT); % number of data batches

fid         = fopen(in_file_path, 'r'); % open for reading raw data
if fid<3
    error('Could not open %s for reading.',ops.fbinary);
end
fidW        = fopen(out_file_path,   'w+'); % open for writing processed data
if fidW<3
    error('Could not open %s for writing.',ops.fproc);    
end
twind=0;
ntb = 64;
w_edge = linspace(0, 1, ntb)';
datr_prev = zeros(ntb, NchanTOT, 'single');

for ibatch = 1:Nbatch
    % we'll create a binary file of batches of NT samples, which overlap consecutively on ops.ntbuff samples
    % in addition to that, we'll read another ops.ntbuff samples from before and after, to have as buffers for filtering
    offset = max(0, twind + 2*NchanTOT*(NT * (ibatch-1) - ntb)); % number of samples to start reading at.
    
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file

    buff = fread(fid, [NchanTOT, NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    nsampcurr = size(buff,2); % how many time samples the current batch has
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); % pad with zeros, if this is the last batch
    end
    if offset==0
        bpad = repmat(buff(:,1), 1, ntb);
        buff = cat(2, bpad, buff(:, 1:NTbuff-ntb)); % The very first batch has no pre-buffer, and has to be treated separately
    end

    datr = single(buff');
    datr = datr - mean(datr,1);
    [spike_row, spike_col, spike_amp] = find_spikes(datr,'threshold',-100);
    datr = clip_out_spikes(datr, spike_row, spike_col, spike_amp);
    
    % datr    = gpufilter(buff, ops, chanMap); % apply filters and median subtraction
    
%     datr(ntb + [1:ntb], :) = datr_prev;
    datr(ntb + [1:ntb], :) = w_edge .* datr(ntb + [1:ntb], :) +...
        (1 - w_edge) .* datr_prev;
   
    datr_prev = datr(ntb +NT + [1:ntb], :);
    datr    = datr(ntb + (1:NT),:); % remove timepoints used as buffers

    datcpu  = int16(datr'); % convert to int16, and gather on the CPU side
    count = fwrite(fidW, datcpu, '*int16'); % write this batch to binary file
    if count~=numel(datcpu)
        error('Error writing batch %g to %s. Check available disk space.',ibatch,ops.fproc);
    end
end
fclose(fidW); % close the files
fclose(fid);

