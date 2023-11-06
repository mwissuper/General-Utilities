function indHS = getHS2(markerPos)

% markerPos is vertical position of marker during period of interest
% (walking forward or backward). Use whole trial of data. Don't look for
% last peak before end of trial

% Get index of vertical markerPos where HS occurred based on extrema. 
x = 1:length(markerPos);
[p,ind] = findpeaks(markerPos(1:end),x,'MinPeakProminence',0.01);
indHS = [];
% Find minima between peaks
for i = 2:length(ind)
    [m(i-1),temp] = min(markerPos(ind(i-1):ind(i)));
    indHS(i-1) = temp + ind(i-1) - 1;
end


% plot(markerPos);hold on,plot(indHS,m,'x'),plot(ind',p,'o');

