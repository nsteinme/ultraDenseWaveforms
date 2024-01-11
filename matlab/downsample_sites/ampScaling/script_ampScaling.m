

% Script to start investigating the effects of downsampling data on
% recorded neurons

% step 1 is to downsample the neuron with different site sizes and
% characterize amplitude and SNR of resulting spikes

%% paths

addpath(genpath(fullfile(githubDir,'ultraDenseWaveforms')))

addpath(genpath(fullfile(githubDir,'npy-matlab')))
%% load data

% dataDir = 'D:\NPUltraWaveforms';
dataDir = '/Users/nicksteinmetz/Dropbox/projects/ultradense/NPUltraWaveforms';

xc = readNPY(fullfile(dataDir, 'channels.xcoords.npy')); 
yc = readNPY(fullfile(dataDir, 'channels.ycoords.npy')); 
wfs = readNPY(fullfile(dataDir, 'clusters.waveforms.npy')); 
ccfCoords = readNPY(fullfile(dataDir, 'clusters.CCF_APDVLR.npy')); 

%%
unitTab = readtable(fullfile(githubDir, 'analysis', 'npUltra_resampleLoss',...
    'allClusterNew3_new-algorithm_python_MB_artefacts.csv')); 

excl = unitTab.artefacts; 

wfs = wfs(~excl,:,:); 
unitTab = unitTab(~excl,:); 

fp = unitTab.footprint; 


%%
sitesize = [12 11 10 9   8   7   6   5   4    3    2    1];
enoise = [4.6 4.9 5.2 5.5 5.9 6.3 6.8 7.6 8.5 9.9 11.5 18];
nADC = 4.5;
tnoise = (enoise.^2 + nADC^2).^(0.5);

upsampRes = 1; % microns


%% 

% cidx = 102; % the large one from the original figure
cidx = 151; % this is a good small one
% cidx = 4498;
thisWF = squeeze(wfs(cidx,:,:))';

wfUp = upsampleWF(thisWF, xc, yc, upsampRes);

% test plot

figure; 
im = imagesc(wfUp(:,:,1)');
caxis(max(abs(wfUp(:)))*[-1 1]);
colormap(colormap_RedWhiteBlue)
axis image; colorbar;
nt = size(wfUp,3); 
for q = 1:nt*10
    im.CData = wfUp(:,:,mod(q-1,nt)+1)'; 
    drawnow; 
    pause(floor(q/10)*0.005);
end

%% test version of upsampling
% this code is assuming that upsampRes is 1 um!


[~,pkX] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],2),[],1);
[~,pkY] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],1),[],2);

mxSz = 10;
plotSamps = 30:70;
tWF = (plotSamps-plotSamps(1))/30;
colors = copper(mxSz+3);

f = figure; f.Color = 'w'; f.Position = [1000        1052         933         286];
subplot(1,3,1);
for ss = 1:mxSz
    idxX = pkX-round(ss/2)+1:(pkX-round(ss/2)+ss-1)+1;
    idxY = pkY-round(ss/2)+1:(pkY-round(ss/2)+ss-1)+1;
    dswf = squeeze(mean(mean(wfUp(idxX, idxY, :), 2), 1));
    wfamp(ss) = max(dswf)-min(dswf);
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
plot((1:mxSz).^2, wfamp./tnoise(end:-1:end-mxSz+1), '.-'); ylabel('SNR'); xlabel('Site size (µm^2)')
box off; 

%% run this on all units

mxSz = 12;
wfamps = nan(size(wfs,1), mxSz); 
snr = nan(size(wfs,1), mxSz); 

for cidx = 1:size(wfs,1)
    if mod(cidx,10)==0
        fprintf(1,'%d...', cidx);
    end
    if mod(cidx,100)==0
        fprintf(1,'\n');
    end
    thisWF = squeeze(wfs(cidx,:,:))';

    wfUp = upsampleWF(thisWF, xc, yc, upsampRes);

    [~,pkX] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],2),[],1);
    [~,pkY] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],1),[],2);
    for ss = 1:mxSz
        idxX = pkX-round(ss/2)+1:(pkX-round(ss/2)+ss-1)+1;
        idxY = pkY-round(ss/2)+1:(pkY-round(ss/2)+ss-1)+1;
        if min(idxX>0) && min(idxY>0) && max(idxX)<=size(wfUp,1) && max(idxY)<=size(wfUp,2)
            dswf = squeeze(mean(mean(wfUp(idxX, idxY, :), 2), 1));
            wfamps(cidx,ss) = max(dswf)-min(dswf);
            
        end
        snr(cidx,:) = wfamps(cidx,:)./tnoise(end:-1:end-mxSz+1);
    end
end

[~,peakSize] = max(snr,[],2);

incl = ~isnan(wfamps(:,mxSz));

%% save calculation results

save bestSiteSizes.mat wfamps snr peakSize


%% plot 

figure; 
subplot(1,3,1); 
hist(peakSize(incl),1:12);

h = histogram2(fp(incl),peakSize(incl),[0:10:100],[4.5:12.5]);

figure; plot(fp(incl), peakSize(incl),'.')
% todo: 
% - assemble final fig

%% final figure

% panels: 
% - zoom of a small WF, overlaid with a few different site sizes
% - example waveforms, amps, SNRs of a small neuron
% - example waveforms, amps, SNRs of a large neuron
% - summary distribution of best size
% - summary distribution of best size as a function of footprint

siteSz = 6; 
mxSz = 12;

colors = myCopper(0.3, mxSz+3);

figure; 


% 1. zoom of small WF


cidx = 2586; % this is a good small one
thisWF = squeeze(wfs(cidx,25:65,:))';
[~,tIdx] = max(max(abs(thisWF)));
wfUp = upsampleWF(thisWF, xc, yc, upsampRes);

subplot(3, 4, [1 5]); 
im = imagesc(wfUp(:,:,tIdx)');
caxis(max(abs(wfUp(:)))*[-1 1]);
colormap(colormap_RedWhiteBlue)
axis image; colorbar;
set(gca, 'YDir', 'normal'); 

hold on; 
ysc = 0.05;
for q = 1:size(thisWF,1)
    plot(xc(q)+linspace(1, siteSz, size(thisWF,2)), ...
        yc(q)+thisWF(q,:)*ysc, 'k', 'LineWidth', 1.5); 
end
title(sprintf('Example %s unit', unitTab.CCFAcronym{cidx}));

 
[~,pkX] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],2),[],1);
[~,pkY] = max(max(max(wfUp,[],3)-min(wfUp,[],3),[],1),[],2);

for ss = [3 7 11]
        idxX = pkX-round(ss/2)+1:(pkX-round(ss/2)+ss-1)+1;
        idxY = pkY-round(ss/2)+1:(pkY-round(ss/2)+ss-1)+1;
        plot([0 0 1 1 0]*ss+idxX(1)-0.5, [0 1 1 0 0]*ss+idxY(1)-0.5, ...
        'Color', colors(ss,:), 'LineWidth', 2.0);
end
xlim([0 30]); ylim([150 250])

%%

figure; 
npUltraPlotWF(xc, yc, thisWF)

%% find the ones with peak wf not on edge sites



chanAmps = squeeze(max(abs(wfs), [], 2)); 
[~, maxCh] = max(chanAmps,[],2); 
maxY = yc(maxCh); maxX = xc(maxCh); 

onProbe = maxY>0 & maxX>0 & maxY<max(yc) & maxX<max(xc);

f = wfbrowse(wfs(onProbe,:,:), xc, yc, unitTab(onProbe,:));


%% try a movie

%280 vpm large
%417 ACB small
onProbeIdx = find(onProbe);

upsampRes = 1; % um
siteSz = 5; 
yRange = [-40 40];
tRange = 25:80;

f = figure; f.Color = 'w';

subplot(1,2,1); cidx = onProbeIdx(280); 

thisWF = squeeze(wfs(cidx,tRange,:))';
[mxAmp,tIdx] = max(max(abs(thisWF)));
[wfUp, upX, upY] = upsampleWF(thisWF, xc, yc, upsampRes);

chanAmps = max(abs(thisWF), [], 2); 
[~, maxCh] = max(chanAmps); 
maxYC = yc(maxCh);

imH = imagesc(upX+3, upY, wfUp(:,:,tIdx)');
axis image;
set(gca, 'YDir', 'normal');
colormap(colormap_RedWhiteBlue)
hold on;
axis off;

caxis(max(abs(wfUp(:)))*[-1 1]/1.5);

% ysc = 0.05;
ysc = 5/mxAmp; 
trH = [];
for q = 1:size(thisWF,1)
    trH(q) = plot(xc(q)+linspace(1, siteSz, size(thisWF,2)), ...
        yc(q)+thisWF(q,:)*ysc, 'k', 'LineWidth', 1.25);
end

ylim(maxYC+yRange);
% h = title(sprintf('%s example neuron', unitTab.CCFAcronym{cidx}));
h = title('Thalamus example neuron'); 
h.FontSize = 14; h.FontName = 'Arial';

% --- 
subplot(1,2,2); cidx = onProbeIdx(417); 

thisWF2 = squeeze(wfs(cidx,tRange,:))';
[mxAmp,tIdx] = max(max(abs(thisWF2)));
[wfUp2, upX, upY] = upsampleWF(thisWF2, xc, yc, upsampRes);

chanAmps = max(abs(thisWF2), [], 2); 
[~, maxCh] = max(chanAmps); 
maxYC = yc(maxCh);

imH2 = imagesc(upX+3, upY, wfUp2(:,:,tIdx)');
axis image;
set(gca, 'YDir', 'normal');
colormap(colormap_RedWhiteBlue)
hold on;
axis off;

caxis(max(abs(wfUp2(:)))*[-1 1]/1.5);

% ysc = 0.05;
ysc = 5/mxAmp; 
trH2 = [];
for q = 1:size(thisWF,1)
    trH2(q) = plot(xc(q)+linspace(1, siteSz, size(thisWF2,2)), ...
        yc(q)+thisWF2(q,:)*ysc, 'k', 'LineWidth', 1.25);
end

ylim(maxYC+yRange);
% h = title(sprintf('%s example neuron', unitTab.CCFAcronym{cidx})); 
h = title('Striatum example neuron'); 
h.FontSize = 14; h.FontName = 'Arial';

vidObj = VideoWriter('exNeurons.mp4', 'MPEG-4');
vidObj.FrameRate =  10;
open(vidObj);

for tt = 1:numel(tRange)
    set(imH, 'CData', wfUp(:,:,tt)');
    set(imH2, 'CData', wfUp2(:,:,tt)');
    drawnow; pause(0.1)
    gf = getframe(f); 
    writeVideo(vidObj,gf);
end
close(vidObj); 