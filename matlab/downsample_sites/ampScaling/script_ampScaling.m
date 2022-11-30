

% Script to start investigating the effects of downsampling data on
% recorded neurons

% step 1 is to downsample the neuron with different site sizes and
% characterize amplitude and SNR of resulting spikes

%% paths

addpath(genpath(fullfile(githubDir,'ultraDenseWaveforms')))
%% load data

dataDir = 'D:\NPUltraWaveforms';

xc = readNPY(fullfile(dataDir, 'channels.xcoords.npy')); 
yc = readNPY(fullfile(dataDir, 'channels.ycoords.npy')); 
wfs = readNPY(fullfile(dataDir, 'clusters.waveforms.npy')); 
ccfCoords = readNPY(fullfile(dataDir, 'clusters.CCF_APDVLR.npy')); 

%%
sitesize = [12 11 10 9   8   7   6   5   4    3    2    1];
enoise = [4.5 5 5.2 5.8 6.0 6.2 6.8 7.8 8.5 10.1 11.5 18];
nADC = 4.5;
tnoise = (enoise.^2 + nADC^2).^(0.5);

upsampRes = 1; % microns


%% 

cidx = 102; 
thisWF = squeeze(wfs(cidx,:,:))';

wfUp = upsampleWF(thisWF, xc, yc, upsampRes);

%% test plot

figure; 
im = imagesc(wfUp(:,:,1)');
caxis(max(abs(wfUp(:)))*[-1 1]);
colormap(colormap_RedWhiteBlue)
axis image; colorbar;
nt = size(wfUp,3); 
for q = 1:nt*10
    im.CData = wfUp(:,:,mod(q-1,nt)+1)'; 
    drawnow; 
end

%% test version of upsampling
% this code is assuming that upsampRes is 1 um!

[~,pkX] = min(min(min(wfUp,[],3),[],2),[],1);
[~,pkY] = min(min(min(wfUp,[],3),[],1),[],2);

mxSz = 10;
plotSamps = 30:70;
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
plot((1:mxSz).^2, wfamp./tnoise(mxSz:-1:1), '.-'); ylabel('SNR'); xlabel('Site size (µm^2)')
box off; 
