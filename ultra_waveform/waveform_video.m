%% Load the results of kilosort. 
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
%Add paths
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
% This uses a function in the github/cortex-lab/spikes repository.
sp = loadKSdir('.') % assumes your current directory has the data in it, if not then replace the directory
%%
sp1 = sp.cids(sp.cgs==2);
xc = (sp.xcoords)/6+1; yc = (sp.ycoords)/6+1;

% unwhiten all the templates
tempsUnW = zeros(size(sp.temps));
for t = 1:size(sp.temps,1)
    tempsUnW(t,:,:) = squeeze(sp.temps(t,:,:))*sp.winv;
end

v = VideoWriter('KSJ_Neighbors96_short.avi');
v.FrameRate=5;

open(v);
t1 = size(tempsUnW,2);


kidx = 21:61; k0 = kidx(1);
imH = [];
for kx = 1:numel(kidx)
    k = kidx(kx);  
    fprintf(['frame' num2str(k-k0) '\n'])
    tic
    %% plot waveform
    % This plot should look like the waveform you see in phy when you press "w" to look at the template. Probably we should use mean waveforms rather than templates though - see function "getWaveforms" in the spikes repository.
    if kx==1
        figure('Position', [10 300 1800 600])
    end
    for i = 1:length(sp1)-1
        if kx==1
            subplot(1,length(sp1),i)
        end
        wf = squeeze(tempsUnW(sp1(1,i)+1,:,:)); % note this was called "62" in phy because of 0-indexing
        nSamp = size(wf, 1); xwf = [1:nSamp]/nSamp;

        %% Tur n the waveform into an image by re-arranging the channels.
        xc = (sp.xcoords)+1; yc = (sp.ycoords)+1;
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
        %
        A = wfIm(:,:,k)';
        % Vq = interp2(A,1);
        % figure('Renderer', 'painters', 'Position', [10 300 300 1200])
        if kx==1
            imH(i) = imagesc(A);

            hold on
            maxx = max((wfIm(:)));
            minx = min((wfIm(:)));
            caxis([minx maxx]); 
            colormap(flipud(colormap_BlueWhiteRed));
        else
            set(imH(i), 'CData', A);
        end
        if kx==1
            
            for ch = 1:numel(sp.xcoords)
                plot(xc(ch)+xwf(20:61)-0.5, yc(ch)+0.5*wf(20:61,ch),'k');
                hold on;
                
            end
            set(gca,'YDir','normal')
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            if i == 1
                title({'Neighbor96',['cluster' num2str(sp1(1,i))]});
            else
                title(['cluster' num2str(sp1(1,i))])
            end
        end
    end
    % fname = ['6micron' num2str(k-k0+1) '.png'];
%     saveas(gcf,fname);
    fname = ['12micron'];
    savefig(fname)
%     frame = getframe(gcf);
%     writeVideo(v,frame);

    close all
    toc
end
% close(v);

