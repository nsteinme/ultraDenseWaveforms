

function newWF = shiftWF(thisWF, ver)

sampShifts = npSampShifts(ver);

assert(numel(sampShifts)==size(thisWF, 1))

fs = 30000;
newWF = zeros(size(thisWF));
for ch = 1:numel(sampShifts)
    newWF(ch,:) = phaseShiftSig(thisWF(ch,:), fs, sampShifts(ch));
end
