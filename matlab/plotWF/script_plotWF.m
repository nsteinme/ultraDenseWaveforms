
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

%%

sitesize = [12 11 10 9   8   7   6   5   4    3    2    1];
enoise = [4.5 5 5.2 5.8 6.0 6.2 6.8 7.8 8.5 10.1 11.5 18];
nADC = 4.5;
tnoise = (enoise.^2 + nADC^2).^(0.5);

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

% inclCID = [cids(1) 26]; 
inclCID = cids;
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
% 8, 15
q = 26; 

%26 is a good small one, so is 30
% 31 is a good big one

makeMovie = false;

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
% absPk = posPk; absPk(negPk<-posPk) = negPk(negPk<-posPk);
absPk = negPk;


% plots

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

nsp = 5;  

f = figure; f.Position = [ 1000         415        1414         923];

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

ax2 = subplot(1, nsp, 2);
cm = parula(200); cm = cm(end:-1:1,:); cm(1,:) = [0.4 0.4 0.4]; 
imagesc(xx(1,:), yy(:,1), negLat');
axis image
colormap(ax2, cm);
ax2.YDir = 'normal';
caxis([min(negLat(:))-1 max(negLat(:))+1]);
axis off;
title('latency to negative peak'); 

ax3 = subplot(1, nsp, 3);
imagesc(xx(1,:), yy(:,1), posLat');
axis image
colormap(ax3, cm);
ax3.YDir = 'normal';
caxis([min(posLat(:))-1 max(posLat(:))+1]);
axis off;
title('latency to positive peak'); 


ax4 = subplot(1,nsp,4); 

im = imagesc(xx(1,:), yy(:,1), wfUp(:,:,90)');
caxis(cax); 
% colormap(ax4, colormap_RedWhiteBlue);
axis image
ax4.YDir = 'normal';
axis off

if ~makeMovie
    t1 = clock;
    updateIm = @(x, dat, tstart)set(x, 'CData', dat(:,:,max(1,ceil(mod(etime(clock, tstart)/5,1)*size(dat,3))))');
    tim = timer('ExecutionMode', 'fixedSpacing', 'Period', 1/30, ...
        'TimerFcn', 'updateIm(im, wfUp(:,:,plotSamps), t1)');
    start(tim);
end


ax5 = subplot(1, nsp, 5); 
im2 = imagesc(reshape(thisWF(:,90), 8, 48)'); 
caxis(cax); 
axis image; 
ax5.YDir = 'normal'; 
axis off; 

if ~makeMovie
    t1 = clock;
    updateIm2 = @(x, dat, tstart)set(x, 'CData', dat(:,:,max(1,ceil(mod(etime(clock, tstart)/5,1)*size(dat,3))))');
    tim2 = timer('ExecutionMode', 'fixedSpacing', 'Period', 1/30, ...
        'TimerFcn', 'updateIm(im2, wfR(:,:,plotSamps), t1)');
    start(tim2);
end

if makeMovie
    outFile = fullfile(pwd, ['wfMovie_' num2str(q)]);
    v = VideoWriter(outFile, 'MPEG-4'); 
    open(v);
    
    for xx = 80:110
        set(im, 'CData', wfUp(:,:,xx)'); 
        set(im2, 'CData', wfR(:,:,xx)'); 
        fr = getframe(f); 
        writeVideo(v,fr);
    end
    close(v);
end


%% test waveform shifting

addpath('C:\Users\nicks\Dropbox\code\analysis\npultra_shiftwf')
nsp = 7;

tic; wfS = shiftWF(thisWF,1); toc;
wfSR = permute(reshape(wfS', size(thisWF,2), 8, 48), [2 3 1]);

plotSamp = 90;

 figure ;
 subplot(1,nsp,1);
 imagesc(wfR(:,:,plotSamp)'); 
 axis image; caxis(cax); axis off; colormap(colormap_RedWhiteBlue)
 set(gca, 'YDir', 'normal');
 
 subplot(1,nsp,2); hold on;
 
 for ch = 1:size(thisWF, 1)
    onewf = thisWF(ch,plotSamps); 
    
    plot(sp.xcoords(ch)+wfT(plotSamps)*tScale, ...
        sp.ycoords(ch)+thisWF(ch,plotSamps)*yScale, ...
        'k'); 
 end
axis image; axis off;
 
subplot(1,nsp,3);
imagesc(wfSR(:,:,plotSamp)'); 
 axis image; caxis(cax); axis off; colormap(colormap_RedWhiteBlue)
 set(gca, 'YDir', 'normal');
 
 subplot(1,nsp,4); hold on;
 
 for ch = 1:size(thisWF, 1)
    onewf = wfS(ch,plotSamps); 
    
    plot(sp.xcoords(ch)+wfT(plotSamps)*tScale, ...
        sp.ycoords(ch)+thisWF(ch,plotSamps)*yScale, ...
        'k'); 
 end
axis image; axis off;

subplot(1, nsp, 5:7); 
xp = 1; yp = 36;
plot(squeeze(wfR(xp,yp,:)), '.-')
hold on;
plot(squeeze(wfSR(xp,yp,:)), '.-')
xlim([plotSamps(1) plotSamps(end)])


%%

idx = 1; plotSamp = plotSamps(idx);

figure; 

subplot(1,2,1);
    im1= imagesc(wfR(:,:,plotSamp)'); 
 axis image; caxis(cax); axis off; colormap(colormap_RedWhiteBlue)
 set(gca, 'YDir', 'normal');

 subplot(1,2,2);
im2 = imagesc(wfSR(:,:,plotSamp)'); 
 axis image; caxis(cax); axis off; colormap(colormap_RedWhiteBlue)
 set(gca, 'YDir', 'normal');
 
while 1
    
    im1.CData = wfR(:,:,plotSamp)';
    im2.CData = wfSR(:,:,plotSamp)';
idx = idx+1; 
if idx>numel(plotSamps); idx = 1; end
plotSamp = plotSamps(idx);

drawnow;
pause(1/7)
end


%% larger and larger recording sites

[~,pkX] = min(min(min(wfUp,[],3),[],2),[],1);
[~,pkY] = min(min(min(wfUp,[],3),[],1),[],2);

mxSz = 12;
plotSamps = 80:120;
tWF = (plotSamps-plotSamps(1))/30;
colors = copper(mxSz+3);

f = figure; f.Color = 'w'; f.Position = [1000        1052         933         286];
subplot(1,3,1);
for ss = 1:mxSz
    idxX = pkX-round(ss/2):(pkX-round(ss/2)+ss-1);
    idxY = pkY-round(ss/2):(pkY-round(ss/2)+ss-1);
    dswf = squeeze(mean(mean(wfUp(idxX, idxY, :), 2), 1));
    wfamp(ss) = -min(dswf);
    plot(tWF, dswf(plotSamps), 'Color', colors(ss,:)); hold on; 
end
xlabel('Time (ms)'); ylabel('Voltage (µV)'); h = legend(array2stringCell((1:mxSz).^2));
h.FontSize = 6;
box off; 

subplot(1,3,2); 
plot((1:mxSz).^2, wfamp, '.-'); ylabel('WF amplitude (µV)'); xlabel('Site size (µm^2)')
ylim([0 max(ylim())]); 
box off; 

subplot(1,3,3); 
plot((1:mxSz).^2, wfamp./tnoise(end:-1:1), '.-'); ylabel('SNR'); xlabel('Site size (µm^2)')
box off; 

print(f, sprintf('exWaveform%d', q), '-dpdf', '-bestfit', '-painters')