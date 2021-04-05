

function sampShifts = npSampShifts(ver)

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

sampShifts = (sampShifts-nChPerADC/2)/nChPerADC;
    