function [posWorkTot, negWorkTot, posWorkMean, negWorkMean] = getPosNegWork(power,fs)

% Calculate positive work as area under power curve where power is positive
% Calculate negative work as area above power curve where power is negative

% MW 11/1/19

indPos = find(power > 0);
posWorkTot = nansum(power(indPos))./fs;
posWorkMean = nanmean(cumsum(power(indPos))./fs);
indNeg = find(power < 0);
negWorkTot = nansum(power(indNeg))./fs;
negWorkMean = nanmean(cumsum(power(indNeg))./fs);
