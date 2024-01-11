

function npUltraPlotWF(xcoords, ycoords, wf)
% wf is channels x time points, in units of uV
% xc and yc are in units of microns

% scale bars are 10um (horiz) and 1000 uV (vert)

Fs = 30000;
upsampRes = 2; % microns
tScale = 2.5e3; yScale = 5;

pk2pk = max(wf,[],2)-min(wf,[],2); 
mxPk2pk = max(pk2pk);
wf = wf/mxPk2pk;

nT = size(wf, 2); 
wfRange = (1:nT)-round(nT/2);
wfT = (wfRange(1):wfRange(end))/Fs;

xc = unique(xcoords); yc = unique(ycoords); 
upX = xc(1):upsampRes:xc(end);
upY = yc(1):upsampRes:yc(end);
[xx,yy] = meshgrid(upX, upY);

% upsample the waveform at each time point
wfR = permute(reshape(wf', nT, 8, 48), [2 3 1]);
wfUp = zeros(numel(upX), numel(upY), size(wfR, 3)); 
for tt = 1:nT       
    wfUp(:,:,tt) = interp2(xc,yc,wfR(:,:,tt)',xx,yy)';
end

[posPk, posLat] = max(wfUp, [], 3);
% posLat(posPk<max(posPk(:))/5) = NaN;
[negPk, negLat] = min(wfUp, [], 3); 
% negLat(negPk>min(negPk(:))/5) = NaN;
% absPk = posPk; absPk(negPk<-posPk) = negPk(negPk<-posPk);
absPk = negPk;

imagesc(xx(1,:), yy(:,1), absPk');
caxis([-1 1]*max(abs(caxis())));
cax = caxis();
colormap(colormap_RedWhiteBlue);
axis image

hold on; 
for ch = 1:size(wf, 1)
    onewf = wf(ch,:); 
    
    plot(xcoords(ch)+wfT*tScale, ...
        ycoords(ch)+wf(ch,:)*yScale, ...
        'k'); 
end
axis off; 
set(gca, 'YDir', 'normal');
plot([-5 5], [-5 -5], 'k'); % horiz scale bar is space, 10 um

% for a voltage scale bar, the distance between adjacent channels will be
% dist/yScale*mxPk2pk. I.e. if dist between channels is 6, yScale is
% 4, and mxPk2pk is 100 uV, then distance between channels is 150 uV
y100 = 1000/mxPk2pk*yScale; 

plot([-5 -5], [0 1]*y100-5, 'k'); 
