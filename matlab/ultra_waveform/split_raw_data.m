githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%% Add paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\Kilosort2_J')) % path to kilosort folder
%%
% find the binary file
mn = 'ZYE_0007';
td = '2020-09-23';
en = 1;
imec = 0;
probeName = 'p1';
serverRoot = expPath(mn, td, en);

%%
ops.fbinary = getProbeFile(serverRoot, probeName, imec);

ops.NchanTOT    = 385; % total number of channels in your recording
ops.fs = 30000; 
ops.ts = 600; %duration of epoch in seconds
ops.trange = [0 Inf]; % time range to sort
%% 
for ibatch = 1:6
    t0 = 0.2*10^7;
    t3 = 2.3*10^7; %interval time between probe move
    tstart  = round((t0+t3*(ibatch-1))/ops.fs);
    ops.fproc = ['batch' num2str(ibatch) '.dat'];
    tic
    NT       = ops.fs*ops.ts; % number of timepoints per batch
    NchanTOT = ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc
    bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
    nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
    ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
    ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
    ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
    ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start

    fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);

    fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
    fidW        = fopen(ops.fproc,   'w'); % open for writing processed data

    % offset = max(0, 2*NchanTOT*(NT * (ibatch-1))); % number of samples to start reading at.
    offset = 2*NchanTOT*30000*tstart; % number of samples to start reading at.
    % offset = 2*NchanTOT*30000*10; % number of samples to start reading at.
    % if offset==0
    %     ioffset = 0; % The very first batch has no pre-buffer, and has to be treated separately
    % else
    %     ioffset = ops.ntbuff;
    % end
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file

    buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)

    fwrite(fidW, buff, 'int16'); % write this batch to binary file

    fclose(fidW); % close the files
    fclose(fid);

    fprintf('Time %3.0fs. Finished preprocessing %d batches. \n', toc, ibatch);
end