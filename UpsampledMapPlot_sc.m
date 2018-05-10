clearvars; clc

%% Load the results of kilosort.
% This uses a function in the github/cortex-lab/spikes repository.
sp = loadKSdir('.'); % assumes your current directory has the data in it, if not then replace the directory

muaidx = find(sp.cgs == 1);

clu_idx = sp.cids + 1; % phy 0-indexing + 1 -> matlab-indexing
clu_idx(muaidx) = [];

% To get value in µV: temp_unfiltered_uV = temp_unfiltered * scale_uV * voltage_step_size;
voltage_step_size = 0.195e-6;
scale_uV = 1000000;

% bad_channels = [7, 8, 19, 28, 30, 32, 34, 94, 131, 137, 148, 149, 151, 162, 171, 232] + 1;

load('wf_sorted.mat'); % Run 'GetWFs.m' first

HPcutoff = 300;
[b,a] = butter(3,HPcutoff/(sp.sample_rate/2),'high');

xc = sp.xcoords; yc = sp.ycoords;

Nrow = max(xc); % = 17
Ncol = max(yc); % = 15
Nsites = Nrow * Ncol; % = 255

% site 32 is reference 
xc(32) = []; yc(32) = [];

UpSampledrow = length(6:1:Nrow*6);
UpSampledcol = length(6:1:Ncol*6);

[X,Y] = meshgrid((6:6:Ncol*6),(6:6:Nrow*6));
[Xq,Yq] = meshgrid((6:1:Ncol*6),(6:1:Nrow*6)); 
    
nClu = length(clu_idx);
AreaOver50 = zeros(nClu,1);

amp_threshold = 50; % 50µV
waveform255 = struct('wf',[],'amplitude',[]); % before upsampling
InterpWF = struct('interpolatedwf',[],'spkamplitude',[]);

% newcolormap = colorcet('R2');
newcolormap = colorcet('I3'); % default: 256 bins

for ww = 1 : nClu
% for ww = 1
    
    waveform = wf.waveFormsMean(ww,:,:);
    waveform = squeeze(waveform);
    waveform = waveform' * scale_uV * voltage_step_size;
    waveform_f = filtfilt(b,a,waveform);

    % interpolate bad sites by neighbouring good sites
    waveform_f(:,20) = mean([waveform_f(:,16),waveform_f(:,24),waveform_f(:,75),waveform_f(:,102),waveform_f(:,97)],2);
    waveform_f(:,35) = mean([waveform_f(:,227),waveform_f(:,131),waveform_f(:,125),waveform_f(:,65),waveform_f(:,63)],2);
    waveform_f(:,150) = mean([waveform_f(:,209),waveform_f(:,244),waveform_f(:,251),waveform_f(:,179),...
        waveform_f(:,145),waveform_f(:,156),waveform_f(:,151),waveform_f(:,185)],2);
    waveform_f(:,172) = mean([waveform_f(:,171),waveform_f(:,177),waveform_f(:,178),waveform_f(:,175),waveform_f(:,156)],2);
    waveform_f(:,95) = mean([waveform_f(:,40),waveform_f(:,22),waveform_f(:,23),waveform_f(:,45),waveform_f(:,39)],2);
    waveform_f(:,9) = mean([waveform_f(:,2),waveform_f(:,3),waveform_f(:,6),waveform_f(:,13),waveform_f(:,17)],2);
    waveform_f(:,8) = mean([waveform_f(:,60),waveform_f(:,54),waveform_f(:,15),waveform_f(:,58),waveform_f(:,11),...
        waveform_f(:,1),waveform_f(:,5),waveform_f(:,12)],2);
    waveform_f(:,233) = mean([waveform_f(:,164),waveform_f(:,213),waveform_f(:,237),waveform_f(:,207)],2);
    waveform_f(:,163) = mean([waveform_f(:,164),waveform_f(:,213),waveform_f(:,212),waveform_f(:,218)],2);
    waveform_f(:,138) = mean([waveform_f(:,144),waveform_f(:,190),waveform_f(:,135),waveform_f(:,143),waveform_f(:,167),...
        waveform_f(:,187),waveform_f(:,139),waveform_f(:,134)],2);
    waveform_f(:,152) = mean([waveform_f(:,158),waveform_f(:,180),waveform_f(:,157),waveform_f(:,184),waveform_f(:,177),...
        waveform_f(:,153),waveform_f(:,148)],2);
    waveform_f(:,149) = mean([waveform_f(:,154),waveform_f(:,183),waveform_f(:,185),waveform_f(:,180),waveform_f(:,144),...
        waveform_f(:,152),waveform_f(:,184),waveform_f(:,143)],2);
    waveform_f(:,182) = mean([waveform_f(:,240),waveform_f(:,249),waveform_f(:,254),waveform_f(:,155),waveform_f(:,146),...
        waveform_f(:,154),waveform_f(:,183),waveform_f(:,186)],2);
    waveform_f(:,29) = mean([waveform_f(:,224),waveform_f(:,38),waveform_f(:,59),waveform_f(:,225),waveform_f(:,25),...
        waveform_f(:,64),waveform_f(:,60)],2);
    waveform_f(:,132) = mean([waveform_f(:,222),waveform_f(:,27),waveform_f(:,192),waveform_f(:,126),waveform_f(:,168),...
        waveform_f(:,165),waveform_f(:,91)],2);
    waveform_f(:,31) = mean([waveform_f(:,220),waveform_f(:,28),waveform_f(:,222),waveform_f(:,27),waveform_f(:,192),...
        waveform_f(:,132),waveform_f(:,126)],2);
    waveform_f(:,33) = mean([waveform_f(:,219),waveform_f(:,225),waveform_f(:,29),waveform_f(:,228),waveform_f(:,64),...
        waveform_f(:,220),waveform_f(:,28)],2);
    waveform_f(:,34) = mean([waveform_f(:,228),waveform_f(:,33),waveform_f(:,64),waveform_f(:,220),waveform_f(:,28),...
        waveform_f(:,222),waveform_f(:,31),waveform_f(:,27)],2);

    waveform_f(:,32) = []; 

    nSamp = size(waveform_f,1); % = 82

    tointerp = zeros(Nrow,Ncol,nSamp);
    interpolatedwf = zeros(UpSampledrow,UpSampledcol,nSamp);

    for d = 1 : nSamp

        for i = 1:length(xc)

            tointerp(xc(i),yc(i),d) = waveform_f(d,i);

        end

        interpolatedwf(:,:,d) = interp2(X,Y,tointerp(:,:,d),Xq,Yq);

    end
    
    waveform255(ww).wf = tointerp;
    InterpWF(ww).interpolatedwf = interpolatedwf;

    amp_3d = max(interpolatedwf, [], 3) - min(interpolatedwf, [], 3);
    
    InterpWF(ww).spkamplitude = amp_3d;
    
    waveform255(ww).amplitude = max(tointerp, [], 3) - min(tointerp, [], 3);

    [idx, idy] = find(amp_3d >= amp_threshold);
    [~,I] = max(amp_3d(:));
    [I_row, I_col] = ind2sub(size(amp_3d),I);
    
    [Map_outline,~] = boundary(idy,idx,0);

    AreaOver50(ww) = numel(idx);
    
    wfamp = max(waveform_f,[],1) - min(waveform_f,[],1);

    scaledwf = zeros(nSamp,Nsites);

    constantfactor = max(wfamp)/9; % maximum channel scaled to 90% of its amplitude

    for i = 1 : length(wfamp)
        scaledwf(:,i) = waveform_f(:,i)/(wfamp(i)+constantfactor);
    end
    
    coloridx = zeros(1,Nsites);
    scaling = max(scaledwf,[],1) - min(scaledwf,[],1);
    binsz = (max(scaling) - min(scaling))/(size(newcolormap,1)-1);
    
    for j = 1 : length(wfamp)
        coloridx(j) = round((scaling(j) - min(scaling))/binsz)+1;
    end

    xwf = (1:nSamp)/100;
    
    ycc = max(sp.xcoords) - sp.xcoords + 1; ycc(32) = []; xcc = sp.ycoords; xcc(32) = [];
    
    figure;
    hFig = gcf;
    set(hFig,'Position', [50 100 1400 380]); set(hFig,'color','w'); set(hFig, 'Visible', 'off');

    ax3 = subplot(1,3,3);
    set(ax3,'Position',[0.7 0.15 0.25 0.75]);
    imagesc(amp_3d);
    colormap(ax3,parula);  
    hcb3 = colorbar; ylabel(hcb3,'µV'); set(hcb3,'FontSize',12); 
    hold on;plot(I_col,I_row,'kx');
    plot(idy(Map_outline),idx(Map_outline),'k','Linewidth',1);
    xlabel('µm');ylabel('µm');title(['Area >= 50µV: ',num2str(AreaOver50(ww)),' (µm^2)']);
    set(gca,'FontSize',14)
    
    ax2 = subplot(1,3,2); 
    set(ax2,'Position',[0.38 0.15 0.25 0.75]);
    box on; hold on;
    for ch = 1:size(scaledwf,2)
        
        plot(xcc(ch)+xwf-0.4, ycc(ch)+scaledwf(:,ch),...
            'Color',newcolormap(coloridx(ch),:),'Linewidth',1); 
        
    end
    colormap(ax2,newcolormap); hcb2 = colorbar; ylabel(hcb2,'scaling'); set(hcb2,'XTick',0:0.25:1); set(hcb3,'FontSize',12);   
    xlabel('site');ylabel('site');
    set(gca,'FontSize',14)
    set(ax2,'XLim',[0.5 15.5],'YLim',[0.5 17.5]);    
    
    ax1 = subplot(1,3,1); 
    set(ax1,'Position',[0.1 0.15 0.225 0.75]);
    box on; hold on;
    for ch = 1:size(scaledwf,2)
        
        plot(xcc(ch)+xwf-0.4, ycc(ch)+0.9*waveform_f(:,ch)./max(wfamp),'k','Linewidth',1); 
        
    end
    xlabel('site');ylabel('site');
    set(gca,'FontSize',14)
    set(ax1,'XLim',[0.5 15.5],'YLim',[0.5 17.5]);
    
    saveas(hFig,['Single_unit_',num2str(ww)],'png');
    
    ww

end 

save('AreaOver50.mat','AreaOver50');
save('InterpolatedWF.mat','InterpWF');
save('wf255.mat','waveform255');
