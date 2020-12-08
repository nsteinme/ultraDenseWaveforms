%% Load the results of kilosort. 
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%Add paths
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Neuropixel_ultra'));

% This uses a function in the github/cortex-lab/spikes repository.
folder = 'E:\ephys\ZYE_0017\2020-11-19\1\imec0\';
% nsp_max = max(nsp);
%%
for j = 1:16
% subfolder = ['ibatch' num2str(j)];
% sp = loadKSdir([folder subfolder]) % assumes your current directory has the data in it, if not then replace the directory
subfolder  = [folder num2str(j)];
sp = loadKSdir(subfolder);
sp1 = sp.cids(sp.cgs==2);

%% run this, if it is first time
gwfparams.dataDir = subfolder;
gwfparams.fileName = sp.dat_path(3:end-1);
gwfparams.dataType = sp.dtype;
gwfparams.nCh = sp.n_channels_dat;
gwfparams.wfWin = [-40 41];
gwfparams.nWf = 100; % randomly choose 300 spikes from each cluster
gwfparams.spikeTimes = sp.st;
gwfparams.spikeClusters = sp.clu;
gwfparams.Fs = sp.sample_rate;
% In wf, wf.waveformsMean: nClu x 256 (sites) x 82 (nSamp)
wf1 = getWaveForms_ultra(gwfparams);
% save('wf.mat','wf1','-v7.3');
save([subfolder '\wf_sorted.mat'],'wf1','-v7.3');
%%
load([subfolder '\wf_sorted.mat']);
xc = (sp.xcoords)/6+1; yc = (sp.ycoords)/6+1;

% unwhiten all the templates
tempsUnW = zeros(size(sp.temps));
for t = 1:size(sp.temps,1)
    tempsUnW(t,:,:) = squeeze(sp.temps(t,:,:))*sp.winv;
end

t1 = size(tempsUnW,2);

kidx = 41;
k = kidx;  
%% plot waveform
% This plot should look like the waveform you see in phy when you press "w" to look at the template. Probably we should use mean waveforms rather than templates though - see function "getWaveforms" in the spikes repository.
nsp_max = length(sp1);
a = ceil(nsp_max./10);
b = 10;
[a1,b1] = ismember(wf1.unitIDs,sp1);
indx = find(a1);
wf2 = wf1.waveFormsMean(indx,:,:);
wf2 = permute(wf2,[1,3,2]);
wft2 = wf1.spikeTimeKeeps(indx,:);
%%
figure('Position', [10 300 900 a*300])
%%
for i = 1:length(sp1)   
    subplot(a,b,i)
    %%
    wf = squeeze(wf2(i,:,:));
    wft = squeeze(wft2(i,:));
    nwft = sum(~isnan(wft));
    %% reduce waveform baseline offset
    wfBaseline = mean(wf(1:10,:));
    wfBaseline2 = repmat(wfBaseline,82,1);
    wf = wf-wfBaseline2;
    %%
    wfmin = abs(min(wf(:)));
    wfmax2 = abs(max(wf(:)));
    wfmax = max(wfmin,wfmax2);
%     wf = wf/wfmax;
    %%
    nSamp = size(wf, 1); xwf = [1:nSamp]/nSamp;
    %% Turn the waveform into an image by re-arranging the channels.
%     xc = (sp.xcoords)+1; yc = (sp.ycoords)+1;
%         for ch = 1:numel(sp.xcoords)
%             plot(xc(ch)+xwf, yc(ch)+0.2*wf(:,ch)); 
%             hold on;
%         end
%
    wfIm = zeros(max(xc), max(yc), size(wf,1));
    for ch = 1:numel(xc)
        if xc(ch)>0 && yc(ch)>0
            wfIm(xc(ch), yc(ch),:) = wf(:,ch);
        end
    end
    %%
    A = wfIm(:,:,k)';
    imH(i) = imagesc(A);
    axis image
    hold on
    maxx = max((wfIm(:)));
    minx = min((wfIm(:)));
    caxis([minx maxx]); 
    colormap(flipud(colormap_BlueWhiteRed));
    for ch = 1:numel(sp.xcoords)
        plot(xc(ch)+xwf(21:61)-0.5, yc(ch)+0.5*wf(21:61,ch)*0.05,'k');
        hold on;
    end
    set(gca,'YDir','normal')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    titlename1 = ['cluster' num2str(sp1(1,i))];
    titlename2 = ['amp: ' num2str(round(wfmax)) ' uV'];
    titlename3 = ['nspikes: ' num2str(nwft)];
    title({titlename1,titlename2,titlename3});
end
figname = ['meanwav_batch_' num2str(j)];
savefig(figname)
end