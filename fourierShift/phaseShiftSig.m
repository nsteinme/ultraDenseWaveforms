

function sig1 = phaseShiftSig(sig, fs, nSamples)
% function sig1 = phaseShiftSig(sig, fs, nSamples)
%
% Shift the fourier phase components of a vector to perform a sub-sample
% shifting of the data, return the data in the time domain.

n = numel(sig);

f = (-n/2:n/2-1)*fs/n; 

% take fft
y = fftshift(fft(sig))/n;

% shift the phase of each sample in a frequency-dependent manner so the
% absolute time shift is constant across frequencies 
y1 = y.*exp(-2*pi*1i*f*nSamples/fs);

% ifft back to time domain
sig1 = real(n*(ifft(ifftshift(y1))));

end