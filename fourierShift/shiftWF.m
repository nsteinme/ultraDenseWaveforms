

function newWF = shiftWF(thisWF, ver)
% thisWF is a matrix of raw neuropixels data nChannels x nSamples
% ver is either 1 or 2, for a 1.0 or 2.0 Neuropixels probe
% newWF has the same size as thisWF but appropriately shifted


% determine how much each channel should be shifted by for this probe
sampShifts = npSampShifts(ver);

assert(numel(sampShifts)==size(thisWF, 1))

% go through each channel and shift that channel alone
fs = 30000;
newWF = zeros(size(thisWF));
for ch = 1:numel(sampShifts)
    newWF(ch,:) = phaseShiftSig(thisWF(ch,:), fs, sampShifts(ch));
end
