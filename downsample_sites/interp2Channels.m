function [buffTemp] = interp2Channels(buff)
%create test data
% buff = randn(384,30000);
% Turn the waveform into an image by re-arranging the channels.
buff = double(buff);
buff2 = reshape(buff,8,48,[]);
clear buff
%% 2D interpolation 
% 1 second of 384 x 30000 samples took 5 seconds to upsample to 1 micron resolution. Any method faster?
% data channels increase from 8x48 to 36x236 channels, which is ~22 times;
tic
[x1,y1] = meshgrid(1:48,1:8);
[xq,yq] = meshgrid(1:0.2:48,1:0.2:8);
% buff3_grid = cellfun(@(x) griddata(x1,y1,squeeze(x),xq,yq),num2cell(buff2,[1 2]),'UniformOutput',false);
buff3_interp = cellfun(@(x) interp2(x1,y1,squeeze(x),xq,yq,'linear'),num2cell(buff2,[1 2]),'UniformOutput',false);
clear buff2
buffTemp = cat(3,buff3_interp{:});
clear buff3_interp
buffTemp = int16(reshape(buffTemp,[],size(buffTemp,3)));

toc
end

%% downsample back to raw data buff, sanity check 

% buff4_interp = reshape(buff4_interp,36,236,[]);
% buffNew_interp = buff4_interp(1:5:end,1:5:end,:);
% % buffDiff = buffNew_interp(:)-buff(:);
% % buffDiff2 = sum(buffDiff(:));
% 
% figure; 
% %original random time point (10)
% % subplot(3,1,1); imagesc(buff2(:,:,10))
% %upsampled time point
% subplot(3,1,2); imagesc(buff4_interp(:,:,10))
% % downsampled time point back to raw config
% subplot(3,1,3); imagesc(buffNew_interp(:,:,10));
% 
% nsites = 36*236;
% fid = fopen('upsample_1micron.dat','r');
% data = fread(fid, [nsites 30000], '*int16');
% fclose(fid);
% data = reshape(data,36,236,[]);
% subplot(3,1,4); imagesc((data(:,:,10)));

