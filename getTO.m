function indTO = getTO(markerPos)

% markerPos is vertical position of marker during period of interest
% (walking forward or backward). 

% Get index of vertical markerPos where HS occurred based on extrema. 

[p,indTO] = findpeaks(-markerPos,1:length(markerPos),'MinPeakProminence',0.01);



