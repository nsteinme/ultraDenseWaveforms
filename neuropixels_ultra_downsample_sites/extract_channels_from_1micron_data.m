%% extract channels from 1micron data

% load raw data
NchanTOT = 36*236;
fpath = 'E:\ephys\ZYE_0006\2020-08-19\3\imec1';
fname = fullfile(fpath,'upsampled_to_1micron.dat');

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
fidW = fopen('downsampled_to_6micron.dat','w'); % open for writing processed data
tic
for ibatch = 1:Nbatch
    fprintf(['processing batch ' num2str(ibatch) '/' num2str(Nbatch) '\n']);
    offset = 2*NchanTOT*NT*(ibatch-1); % number of samples to start reading at.
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file
    buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    buffTemp = reshape(buff,36,236,[]);
    buffTemp2 = buffTemp(1:5:end,1:5:end,:);
    fwrite(fidW, buffTemp2, 'int16');
end   
fclose(fid);
fclose(fidW);

toc

