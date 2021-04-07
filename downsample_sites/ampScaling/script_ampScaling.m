

% Script to start investigating the effects of downsampling data on
% recorded neurons

% step 1 is to downsample the neuron with different site sizes and
% characterize amplitude and SNR of resulting spikes

%%
sitesize = [12 11 10 9   8   7   6   5   4    3    2    1];
enoise = [4.5 5 5.2 5.8 6.0 6.2 6.8 7.8 8.5 10.1 11.5 18];
nADC = 4.5;
tnoise = (enoise.^2 + nADC^2).^(0.5);


