function [indLHS_new, indRHS_new] = remRepHS(indLHS, indRHS)

% Check if two consec HS events on one side. If yes, remove
% first one.
indLHS_new = indLHS;
indRHS_new = indRHS;
indHS(:,1) = [indLHS'; indRHS'];
indHS(:,2) = [ones(length(indLHS),1); 2*ones(length(indRHS),1)];
[y,iSort] = sort(indHS(:,1));
sideSort = indHS(iSort,2);
sideOff = [sideSort(2:end); nan];
irep = find(sideSort == sideOff);
for i = 1:length(irep)
    rem = find(indLHS_new == y(irep(i)),1,'first');
    indLHS_new(rem) = [];
    rem = find(indRHS_new == y(irep(i)),1,'first');
    indRHS_new(rem) = [];
end