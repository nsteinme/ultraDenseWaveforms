
function [wfUp, upX, upY] = upsampleWF(thisWF, xcoords, ycoords, upsampRes)
% function wfUp = upsampleWF(thisWF, xcoords, ycoords, upsampRes)
% thisWF is chans x samples
xc = unique(xcoords); yc = unique(ycoords); 
upX = xc(1):upsampRes:xc(end);
upY = yc(1):upsampRes:yc(end);
[xx,yy] = meshgrid(upX, upY);

% upsample the waveform at each time point
wfR = permute(reshape(thisWF', size(thisWF,2), numel(xc), numel(yc)), [2 3 1]);
wfUp = zeros(numel(upX), numel(upY), size(wfR, 3)); 
for tt = 1:size(thisWF,2)        
    wfUp(:,:,tt) = interp2(xc,yc,wfR(:,:,tt)',xx,yy,'linear')';
    % tested cubic briefly and it didn't seem to make any difference
end