

%%

fs = 100;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*15*t - pi/4) - sin(2*pi*40*t);


%%

y = fft(x);
% z = fftshift(y);

ly = length(y);
% f = (-ly/2:ly/2-1)/ly*fs;
f = (0:ly-1)/ly*fs;

% stem(f,abs(z))
stem(f,abs(y))
xlabel 'Frequency (Hz)'
ylabel '|y|'
grid
xlim([0 fs/2])

%%

tol = 1e-6;
y(abs(y) < tol) = 0;

theta = angle(y);

stem(f,theta/pi)
xlabel 'Frequency (Hz)'
ylabel 'Phase / \pi'
grid
xlim([0 fs/2])


%% shift and ifft

%%
fs = 1000;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*15*t - pi/4) - sin(2*pi*40*t);
n = numel(t);
delt = 1/fs;
delf = fs/n;
t = (0:n-1)*delt;
f = (-n/2:n/2-1)*delf;
y = fftshift(fft(x))/n;
i0 = n/2+1;
ip = i0+15;
im = i0-15;
first_fund = y(ip)*exp(2*pi*1i*f(ip)*t) + y(im)*exp(2*pi*1i*f(im)*t);
figure;
% plot(t,x,t,first_fund, '.-')              % same as what you have
plot(t,x, '.-') 
% t0 = (angle(y(ip))/(2*pi*f(ip)))

% t0 = pi/32*f(ip);
for t0 = [0:0.02:1]/fs
y1 = y.*exp(-2*pi*1i*f*t0);          % y1(ip) and y1(im) are real.
% y1(1) = real(y1(1));                % band aid fix                
x1 = n*(ifft(ifftshift(y1)));       % approx. circularly shifted signal
first_fund1 = y1(ip)*exp(2*pi*1i*f(ip)*t) + y1(im)*exp(2*pi*1i*f(im)*t);
% figure (5)
hold off;
plot(t,x, '.-')              % same as what you have
hold on;
plot(t,x1, '.-')
drawnow;
end

%%


fs = 1000;
t = 0:1/fs:1-1/fs;
x = cos(2*pi*15*t - pi/4) - sin(2*pi*40*t);

figure; 

for t0 = [0:0.02:1]
    x1 = phaseShiftSig(x,fs,t0);
    hold off; 
    plot(t,x, '.-')             
    hold on;
    plot(t,x1, '.-')
%     xlim([0.25 0.29]);
    drawnow;
end