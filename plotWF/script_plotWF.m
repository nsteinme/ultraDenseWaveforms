
% This is a sandbox script to work out the plotting for NP Ultra waveforms.
% Relevant code moved to functions from here. 

mn = 'ZYE_0007'; 
td = '2020-09-23';
en = 1;
pn = 'imec0';

serverRoot = expPath(mn, td, en);
apFile = getProbeFile(serverRoot, pn);
dataDir = fullfile(apFile);

%% settings


ops = struct(); 
ops.fs = 30000;
ops.fshigh = 300; % high pass filter cutoff
ops.CAR = false;

chanMap = 1:384;

gain = 2.34; % uV/bit

%% 

sp = loadKSdir(dataDir);

incl = sp.st>0; 
st = sp.st(incl);
clu = sp.clu(incl); 

inclCID = [104 108 111 113 117 123]; 

nCID = numel(inclCID);
inclClu = ismember(clu, inclCID); 

%%


gwfparams.dataDir = dataDir;    % KiloSort/Phy output folder
gwfparams.fileName = 'p2_g2_t0.imec0.ap.bin';         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-90 90];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 100;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    round(ops.fs*st(inclClu)); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clu(inclClu); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

wf = getWaveForms(gwfparams);


%% filter and preprocess waveforms
 

%chMeans = mean(mean(wf.waveFormsMean(:,:,1:50),3),1); 
mnSub = zeros(size(wf.waveForms)); 
for u = 1:size(mnSub,1)
    for s = 1:size(mnSub,2)
        %mnSub(u,s,:,:) = squeeze(wf.waveForms(u,s,:,:)) - chMeans';
        buff = squeeze(wf.waveForms(u,s,:,:)); 
        
        if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
            [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
        else
            [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz
        end

        %dataRAW = gpuArray(buff); % move int16 data to GPU
        dataRAW = buff';
        dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
        dataRAW = dataRAW(:, chanMap); % subsample only good channels

        % subtract the mean from each channel
        dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

        % CAR, common average referencing by median
        if getOr(ops, 'CAR', 1)
            dataRAW = dataRAW - median(dataRAW, 2); % subtract median across channels
        end

        % next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
        datr = filter(b1, a1, dataRAW); % causal forward filter
        datr = flipud(datr); % reverse time
        datr = filter(b1, a1, datr); % causal forward filter again
        datr = flipud(datr); % reverse time back

        mnSub(u,s,:,:) = datr';
        
    end
end

%% compute stuff for a neuron

q = 1; 

thisWF = squeeze(mean(mnSub(q,:,:,:),2)); % now it is chans x samples



%% plots


% subplots to make: 
% 1. waveforms with peak map in the background
% 2. latency to negative peak
% 3. latency to positive peak
% 4. movie of waveform

nsp = 4;  

f = figure; 

subplot(1, nsp, 1); 



