sp = loadKSdir('.'); % assumes your current directory has the data in it, if not then replace the directory

gwfparams.dataDir = pwd;
gwfparams.fileName = sp.dat_path;
gwfparams.dataType = sp.dtype;
gwfparams.nCh = sp.n_channels_dat;
gwfparams.wfWin = [-40 41];
gwfparams.nWf = 300; % randomly choose 300 spikes from each cluster
gwfparams.spikeTimes = sp.st;
gwfparams.spikeClusters = sp.clu;
gwfparams.Fs = sp.sample_rate;

% In wf, wf.waveformsMean: nClu x 256 (sites) x 82 (nSamp)
wf = getWaveForms_sc(gwfparams);
% save('wf.mat','wf','-v7.3');
save('wf_sorted.mat','wf','-v7.3');

