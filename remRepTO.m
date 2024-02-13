function [indLTO_new, indRTO_new] = remRepTO(indLTO, indRTO)

% Check if two consec TO events on one side. If yes, remove
% second one.
indLTO_new = indLTO;
indRTO_new = indRTO;
indTO(:,1) = [indLTO'; indRTO'];
indTO(:,2) = [ones(length(indLTO),1); 2*ones(length(indRTO),1)];
[y,iSort] = sort(indTO(:,1));
sideSort = indTO(iSort,2);
sideOff = [sideSort(2:end); nan];
irep = find(sideSort == sideOff) + 1; % get the second one
% Get index values for all repeats
reps = y(irep);
% Remove the repeats from each side's array
[c,i1,i2] = intersect(indLTO_new,reps);
indLTO_new(i1) = [];
[c,i1,i2] = intersect(indRTO_new,reps);
indRTO_new(i1) = [];