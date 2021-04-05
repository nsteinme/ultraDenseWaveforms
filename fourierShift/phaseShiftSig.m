

function sig1 = phaseShiftSig(sig, fs, nSamples)

n = numel(sig);

f = (-n/2:n/2-1)*fs/n;

y = fftshift(fft(sig))/n;

y1 = y.*exp(-2*pi*1i*f*nSamples/fs);

sig1 = real(n*(ifft(ifftshift(y1))));

end