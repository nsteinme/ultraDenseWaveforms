
% This is a sandbox script to work out the plotting for NP Ultra waveforms.
% Relevant code moved to functions from here. 

mn = 'ZYE_0007'; 
td = '2020-09-23';
en = 1;
pn = 'imec0';

% serverRoot = expPath(mn, td, en);
% apFile = getProbeFile(serverRoot, pn);
% dataDir = fullfile(apFile);

apFile = 'D:\DownloadsD\batch2\batch2.dat';
[dataDir, fn, ext] = fileparts(apFile);
fn = [fn ext];

%% settings


ops = struct(); 
ops.fs = 30000;
ops.fshigh = 300; % high pass filter cutoff
ops.CAR = false;

chanMap = 1:384;

gain = 2.34; % uV/bit

upsampRes = 2; % microns

wfRange = [-90 90];
plotSamps = 70:110;
wfT = (wfRange(1):wfRange(end))/ops.fs;

%% 

sp = loadKSdir(dataDir);

incl = sp.st>0; 
st = sp.st(incl);
clu = sp.clu(incl); 

cids = unique(clu); 

inclCID = [cids(1) 26]; 

nCID = numel(inclCID);
inclClu = ismember(clu, inclCID); 

%%


gwfparams.dataDir = dataDir;    % KiloSort/Phy output folder
gwfparams.fileName = fn;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = wfRange;              % Number of samples before and after spiketime to include in waveform
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

        % subtract the mean from each channel, using the "baseline"
        dataRAW = dataRAW - mean(dataRAW(1:20,:), 1); % subtract mean of each channel

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

q = 2; 

xc = unique(sp.xcoords); yc = unique(sp.ycoords); 
upX = xc(1):upsampRes:xc(end);
upY = yc(1):upsampRes:yc(end);
[xx,yy] = meshgrid(upX, upY);


thisWF = gain*squeeze(nanmean(mnSub(q,:,:,:),2)); % now it is chans x samples

% upsample the waveform at each time point
wfR = permute(reshape(thisWF', size(thisWF,2), 8, 48), [2 3 1]);
wfUp = zeros(numel(upX), numel(upY), size(wfR, 3)); 
for tt = 1:size(thisWF,2)        
    wfUp(:,:,tt) = interp2(xc,yc,wfR(:,:,tt)',xx,yy)';
end

[posPk, posLat] = max(wfUp, [], 3);
posLat(posPk<max(posPk(:))/5) = NaN;
[negPk, negLat] = min(wfUp, [], 3); 
negLat(negPk>min(negPk(:))/5) = NaN;
absPk = posPk; absPk(negPk<-posPk) = negPk(negPk<-posPk);


%% plots

tScale = 4e3; yScale = 0.05;

% subplots to make: 
% 1. waveforms with peak map in the background
% 2. latency to negative peak
% 3. latency to positive peak
% 4. movie of waveform
% -- CSD? 
% -- 3D render of space/time waveform? 
% -- slices over time and one spatial dimension? 
% -- contour plot with different shades of blue for positive deflections at
% different timepoints, different shades of red for negative deflections at
% different timepoints

nsp = 2;  

f = figure; 

ax1 = subplot(1, nsp, 1); 


imagesc(xx(1,:), yy(:,1), absPk');
caxis([-1 1]*max(abs(caxis())));
cax = caxis();
colormap(ax1,colormap_RedWhiteBlue);
axis image

hold on; 
for ch = 1:size(thisWF, 1)
    onewf = thisWF(ch,plotSamps); 
    
    plot(sp.xcoords(ch)+wfT(plotSamps)*tScale, ...
        sp.ycoords(ch)+thisWF(ch,plotSamps)*yScale, ...
        'k'); 
end
title(sprintf('pk-pk amplitude = %.1f uV', max(posPk(:)-negPk(:))))

ax1.YDir = 'normal';
axis off

% ax2 = subplot(1, nsp, 2);
% cm = parula(200); cm = cm(end:-1:1,:); cm(1,:) = [0.4 0.4 0.4]; 
% imagesc(xx(1,:), yy(:,1), negLat');
% axis image
% colormap(ax2, cm);
% ax2.YDir = 'normal';
% caxis([min(negLat(:))-1 max(negLat(:))+1]);
% axis off;
% title('latency to negative peak'); 
% 
% ax3 = subplot(1, nsp, 3);
% imagesc(xx(1,:), yy(:,1), posLat');
% axis image
% colormap(ax3, cm);
% ax3.YDir = 'normal';
% caxis([min(posLat(:))-1 max(posLat(:))+1]);
% axis off;
% title('latency to positive peak'); 


ax4 = subplot(1,nsp,2); 

im = imagesc(xx(1,:), yy(:,1), wfUp(:,:,90)');
caxis(cax); 
% colormap(ax4, colormap_RedWhiteBlue);
axis image
ax4.YDir = 'normal';
axis off
t1 = clock; 
updateIm = @(x, dat, tstart)set(x, 'CData', dat(:,:,max(1,ceil(mod(etime(clock, tstart)/5,1)*size(dat,3))))'); 
tim = timer('ExecutionMode', 'fixedSpacing', 'Period', 1/30, ...
    'TimerFcn', 'updateIm(im, wfUp(:,:,plotSamps), t1)'); 
start(tim); 
