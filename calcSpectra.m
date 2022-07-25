function [f, P1] = calcSpectra(x, fs)
   
L = length(x);        % Length of signal
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;