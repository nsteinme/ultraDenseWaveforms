

function sampShifts = npSampShifts(ver)
% return the true sampling time (relative to the average sample time) of
% each channel of a neuropixels 1.0 or 2.0 probe. 

nCh = 384;
sampShifts = zeros(nCh, 1);
switch ver 
    case 1
        nADC = 32; 
        
    case 2
        nADC = 24;
    otherwise
        error('unrecognized probe version');
        return; 
end

nChPerADC = nCh/nADC;
startChan = [1:nChPerADC*2:nCh 2:nChPerADC*2:nCh];
for n = 1:nADC
    sampShifts(startChan(n):2:startChan(n)+2*nChPerADC-1) = 1:nChPerADC;
end

% here just centering the shift amount (arbitrary) and returning in units
% of samples
sampShifts = (sampShifts-nChPerADC/2)/nChPerADC;
    