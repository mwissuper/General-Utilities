function [posPeaks,iPosPeaks,negPeaks,iNegPeaks] = getPosNegPeaks(x)

% Find the extrema for the positive and negative values of x separately.
% Need this for power and force data. Use output to calculate mean of peaks
% and for plotting

% MW 11/1/19

% Positive value peaks
if isempty(find(x>0)) % no positive values
    posPeaks = nan; iPosPeaks = nan;
else
    [posPeaks,iPosPeaks] = findpeaks(x,'MinPeakProminence',0.05);
    % Remove extrema values that are negative
    indRemove = find(posPeaks<0);
    posPeaks(indRemove) = [];
    iPosPeaks(indRemove) = [];
end

% Negative value peaks
if isempty(find(x<0)) % no positive values
    negPeaks = nan; iNegPeaks = nan;
else
    [temp,iNegPeaks] = findpeaks(-x,'MinPeakProminence',0.05);
    negPeaks = -temp;
    % Remove extrema values that are negative
    indRemove = find(negPeaks>0);
    negPeaks(indRemove) = [];
    iNegPeaks(indRemove) = [];
end

% Plot to check peaks are good
% 
% figure
% plot(x),hold on,plot(iPosPeaks,posPeaks,'o'),plot(iNegPeaks,negPeaks,'x')
% xlabel('Sample');



