% buff is matrix of channels x samples
load('neuropixelUltra_kilosortChanMap.mat')

buff = randn(384,3000);
xc = (xcoords)/6+1; yc = (ycoords)/6+1;

%% Turn the waveform into an image by re-arranging the channels.

buff2 = zeros(max(xc), max(yc), size(buff,2));
for ch = 1:numel(xc)
    buff2(xc(ch), yc(ch),:) = buff(ch,:);
end

%% 2D interpolation 
% 1 second of 384 x 30000 samples took 162 seconds to upsample to 1 micron resolution. Any method faster?
% data channels increase from 8x48 to 36x236 channels, which is ~22 times;

% 1 second of 384 x 30000 samples took 50 seconds to upsample to 3 micron resolution. Any method faster?

tic
[x1,y1] = meshgrid(1:48,1:8);
% [xq,yq] = meshgrid(1:0.2:48,1:0.2:8);
[xq,yq] = meshgrid(1:0.5:48,1:0.5:8);
buff3 = cellfun(@(x) griddata(x1,y1,squeeze(x),xq,yq),num2cell(buff2,[1 2]),'UniformOutput',false);
buff4 = cat(3,buff3{:});
toc

%% downsample back to raw data buff, sanity check 
% buffNew = buff4(1:5:end,1:5:end,:);
buffNew = buff4(1:2:end,1:2:end,:);
buffDiff = buffNew(:)-buff(:);
buffDiff2 = sum(buffDiff(:));

figure; 
%original random time point (10)
subplot(3,1,1); imagesc(buff2(:,:,10))
%upsampled time point
subplot(3,1,2); imagesc(buff4(:,:,10))
% downsampled time point back to raw config
subplot(3,1,3); imagesc(buffNew(:,:,10))
