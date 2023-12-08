function indHS = getHS(markerPos)

% markerPos is vertical position of marker during period of interest
% (walking forward or backward). Often bad data at beg of trial, so skip
% that period

% Get index of vertical markerPos where HS occurred based on extrema. 
x = 50:length(markerPos);
[p,ind] = findpeaks(markerPos(50:end),x,'MinPeakProminence',0.01);
indHS = [];
% Find minima between peaks
for i = 2:length(ind)
    [m(i-1),temp] = min(markerPos(ind(i-1):ind(i)));
    indHS(i-1) = temp + ind(i-1) - 1;
end

% For last step, look for min between last peak and end of time
% period of interest. This doesn't work if data is cut off before reach stillness
% Why do I need this last step past the last peak?
[m(i),temp] = min(markerPos(ind(i):end));
indHS(i) = temp + ind(i) - 1;

% plot(markerPos);hold on,plot(indHS,m,'x'),plot(ind',p,'o');

