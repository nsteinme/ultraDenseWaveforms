
%% Load the results of kilosort. 
% This uses a function in the github/cortex-lab/spikes repository.
sp = loadKSdir('.');
sp
% unwhiten all the templates
tempsUnW = zeros(size(sp.temps));
for t = 1:size(sp.temps,1)
    tempsUnW(t,:,:) = squeeze(sp.temps(t,:,:))*sp.winv;
end

%% plot waveform
% This plot should look like the waveform you see in phy when you press "w" to look at the template. Probably we should use mean waveforms rather than templates though - see function "getWaveforms" in the spikes repository.
figure; 
wf = squeeze(tempsUnW(63,:,:)); % note this was called "62" in phy because of 0-indexing
nSamp = size(wf, 1); xwf = [1:nSamp]/nSamp;
xc = sp.xcoords; yc = sp.ycoords;
for ch = 1:numel(sp.xcoords)
    plot(xc(ch)+xwf, yc(ch)+0.2*wf(:,ch)); 
    hold on;
end

%% Turn the waveform into an image by re-arranging the channels.
wfIm = zeros(max(xc), max(yc), size(wf,1));
for ch = 1:numel(xc)
    if xc(ch)>0 && yc(ch)>0
        wfIm(xc(ch), yc(ch),:) = wf(:,ch);
    end
end

%% Play it as a movie.
figure; 
im = imagesc(wfIm(:,:,1)); 
mx = max(abs(wfIm(:))); 
caxis([-mx mx]*0.75); colormap(colormap_BlueWhiteRed);
for s = 1:nSamp
    set(im, 'CData', wfIm(:,:,s)); drawnow; pause(0.05);
end