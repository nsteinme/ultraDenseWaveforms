% sp = loadKSdir('.'); % assumes your current directory has the data in it, if not then replace the directory

gwfparams.dataDir = 'E:\ephys\ZYE_0017\2020-11-19\1\imec0\1';
gwfparams.fileName = sp.dat_path(3:end-1);
gwfparams.dataType = sp.dtype;
gwfparams.nCh = sp.n_channels_dat;
gwfparams.wfWin = [-40 41];
gwfparams.nWf = 100; % randomly choose 300 spikes from each cluster
gwfparams.spikeTimes = sp.st;
gwfparams.spikeClusters = sp.clu;
gwfparams.Fs = sp.sample_rate;
%%
% In wf, wf.waveformsMean: nClu x 256 (sites) x 82 (nSamp)
wf1 = getWaveForms_ultra(gwfparams);
% save('wf.mat','wf1','-v7.3');
save('wf_sorted.mat','wf1','-v7.3');
%%
