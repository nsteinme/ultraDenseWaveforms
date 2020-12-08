githubDir = 'C:\Users\Steinmetz lab\Documents\MATLAB';
%Add paths
addpath(genpath(fullfile(githubDir, 'Neuropixel_ultra')))
%% get binary raw data in workspace
% load filtered data
% NchanTOT = 384;
% fpath = 'E:\ephys\ZYE_0006\2020-08-19\3\imec1\KSJ_Neighbors96';
% fname = fullfile(fpath,'temp_wh.dat');

% load raw data
NchanTOT = 385;
fpath = 'E:\ephys\ZYE_0006\2020-08-19\3\imec1';
fname = fullfile(fpath,'batch1.dat');

% size in bytes of raw binary
cmd = sprintf('stat -Lc %%s %s', fname);
[status, r] = system(cmd);
bytes = str2double(r);
if isnan(bytes)
    o = dir(fname);
    bytes = o.bytes;
end

NT1 = floor(bytes/NchanTOT/2); % number of total timepoints

NT = 64*1024;
Nbatch  = ceil(NT1 /NT); % number of data batches
%%
fid = fopen(fname, 'r'); % open for reading raw data
fidW = fopen('upsampled_to_1micron.dat','w'); % open for writing processed data

for ibatch = 1:Nbatch
    fprintf(['processing batch ' num2str(ibatch) '/' num2str(Nbatch) '\n']); 
    offset = 2*NchanTOT*NT*(ibatch-1); % number of samples to start reading at.
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file
    buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    % channel 385 is sync channel
    buff = buff(1:384,:);
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    [buffTemp] = interp2Channels(buff);
    fwrite(fidW, buffTemp, 'int16');
end   
fclose(fid);
fclose(fidW);



