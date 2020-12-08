%% default options are in parenthesis after the comment

addpath(genpath('C:\Users\Nick\Documents\github\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Nick\Documents\github\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Nick\Documents\github\ultraDenseWaveforms'))

pathToYourConfigFile = 'C:\Users\Nick\Documents\github\ultraDenseWaveforms\kilosortConfig'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'ultraDenseConfig.m'))

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% get list of all filenames
ultraDenseDataDir = '\\zubjects.cortexlab.net\Subjects\kampff\ultra-dense-survey\Data\bin files';
ultraDenseRezDir = '\\zubjects.cortexlab.net\Subjects\kampff\ultra-dense-survey\Data\KSresults'; 
mkdir(ultraDenseRezDir);
fn = dir(fullfile(ultraDenseDataDir, '*.bin'));

for fx = 6:numel(fn)
    try
        fprintf(1, '%d of %d\n', fx, numel(fn));
        
        % set ops to work for this filename
        ops.fbinary = fullfile(ultraDenseDataDir, fn(fx).name);
        saveDir = fullfile(ultraDenseRezDir, fn(fx).name(1:end-4)); mkdir(saveDir);
        
        % run kilosort
        [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
        rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
        rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
        
        % save matlab results file
        save(fullfile(saveDir,  'rez.mat'), 'rez', '-v7.3');
        
        % save python results file for Phy
        rezToPhy(rez, saveDir);
        
        % remove temporary file
        delete(ops.fproc);
    catch me
        disp(me)
        delete(ops.fproc);
    end
end

