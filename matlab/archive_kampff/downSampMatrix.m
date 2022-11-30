

function [downSampGridAlign, downSampGridStag] = downSampMatrix(yOffset, xOffset, makePlots)
% function [downSampGridAlign, downSampGridStag] = downSampMatrix(yOffset, xOffset, makePlots)
%
% Example: 
% >> yOffset = 5; 
% >> xOffset = 5; 
% >> makePlots = true;
% >> [dsAlign, dsStag] = downSampMatrix(yOffset, xOffset, makePlots);


nSitesIncl = 12;
siteSp = 6; 
siteGrid = [17 15];

upSampSize = siteGrid*siteSp; 

np2siteSize = [12 12]; 

np2originsYstag = repmat([0 32 16 48], 1, 100)+yOffset;
np2originsYalign = repmat([0 32], 1, 100)+yOffset;
np2originsX = reshape(repmat(0:15:1000,2,1), 1, [])+xOffset;

np2originsYstag = np2originsYstag(1:nSitesIncl);
np2originsYalign = np2originsYalign(1:nSitesIncl);
np2originsX = np2originsX(1:nSitesIncl);

downSampGridStag = zeros(upSampSize); 

for q = 1:numel(np2originsX)
    downSampGridStag(np2originsX(q)+[0:np2siteSize(1)-1]+1, ...
        np2originsYstag(q)+[0:np2siteSize(2)-1]+1) = q;
end

downSampGridAlign = zeros(upSampSize); 

for q = 1:numel(np2originsX)
    downSampGridAlign(np2originsX(q)+[0:np2siteSize(1)-1]+1, ...
        np2originsYalign(q)+[0:np2siteSize(2)-1]+1) = q;
end

if makePlots
figure; 
subplot(1,2,1);
imagesc(downSampGridStag)
axis image;
set(gca, 'YDir', 'normal');

subplot(1,2,2);
imagesc(downSampGridAlign)
axis image;
set(gca, 'YDir', 'normal');


end