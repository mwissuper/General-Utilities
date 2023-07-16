function ind = findMarkerIndCsv(markerNames,s)

% Find csv data index corresponding to marker names 's' using an array of
% markerNames

for j = 1:length(markerNames)
    idx(j) = strcmp(markerNames(j,:),s);
end
temp = find(idx == 1,1,'first');  
if ~isempty(temp) 
    ind = 2+(temp-1)*3+1; % 1st 2 cols of data are not marker data, then x, y, z components of marker
else
    ind = nan;
    msg = sprintf('%s not present',s);
    disp(msg);
end