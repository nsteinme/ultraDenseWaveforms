
% plan: make it animate by going through timepoints and moving a green dot
% along the traces while updating the image in background
% - save out a series of selected neurons played in sequence

function f = wfbrowse(wfs, xc, yc, unitTab)


% cidx = 2586; % this is a good small one
% cidx = 2705;
cidx = 1;
f = figure; f.Color = 'w';

ud.cidx = cidx;
ud.unitTab = unitTab;
ud.wfs = wfs;
ud.xc = xc;
ud.yc = yc;
ud.imH = [];
ud.trH = []; 
f.UserData = ud;

wfbrowseUpdate(f);

set(f, 'KeyPressFcn', @(f,k)wfbrowseCallback(f, k));


function wfbrowseCallback(f, keydata)

ud = get(f, 'UserData');

switch keydata.Key
    case 'rightarrow'
        ud.cidx = ud.cidx+1;
        if ud.cidx>size(ud.wfs,1)
            ud.cidx = 1;
        end
    case 'leftarrow'
        ud.cidx = ud.cidx-1;
        if ud.cidx<1
            ud.cidx = size(ud.wfs,1);
        end
end

f.UserData = ud;
wfbrowseUpdate(f);

function wfbrowseUpdate(f)
ud = f.UserData; wfs = ud.wfs; xc = ud.xc; yc = ud.yc; cidx = ud.cidx;

upsampRes = 1; % um
siteSz = 5; 

thisWF = squeeze(wfs(cidx,25:80,:))';
[mxAmp,tIdx] = max(max(abs(thisWF)));
[wfUp, upX, upY] = upsampleWF(thisWF, xc, yc, upsampRes);

chanAmps = max(abs(thisWF), [], 2); 
[~, maxCh] = max(chanAmps); 
maxYC = yc(maxCh);

if isempty(ud.imH) 
    ud.imH = imagesc(upX+3, upY, wfUp(:,:,tIdx)');
    axis image; colorbar;
    set(gca, 'YDir', 'normal');
    colormap(colormap_RedWhiteBlue)
    hold on; 
    axis off; 
else
    set(ud.imH, 'CData', wfUp(:,:,tIdx)');
end
caxis(max(abs(wfUp(:)))*[-1 1]);

% ysc = 0.05;
ysc = 5/mxAmp; 
if isempty(ud.trH)
    trH = [];
    for q = 1:size(thisWF,1)
        trH(q) = plot(xc(q)+linspace(1, siteSz, size(thisWF,2)), ...
            yc(q)+thisWF(q,:)*ysc, 'k', 'LineWidth', 1.5); 
    end
    ud.trH = trH;
else
    for q = 1:size(thisWF,1)
        set(ud.trH(q), 'YData', yc(q)+thisWF(q,:)*ysc); 
    end
end
ylim(maxYC+[-40 40]);
title(sprintf('Neuron %d, %s', cidx, ud.unitTab.CCFAcronym{cidx}));

f.UserData = ud;