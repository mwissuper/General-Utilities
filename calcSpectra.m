function [f, P1, P1log, BW] = calcSpectra(x, fs)

% One-sided spectrum up to fNyquist. Also calculate 3dB BW.
   
L = length(x);        % Length of signal
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:floor(L/2))/L; % hz

P1log = 20*log10(P1);
P1log_off = [P1log(2:end); nan];
ind = find(P1log > -3 & P1log_off < -3,1,'first');
BW = f(ind);