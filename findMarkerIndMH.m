function indMarker = findMarkerIndMH(s,markerID)

% Take input string s of name of marker and find which index it is in
% MarkerID. Check for cases of 1 or 2 subjects

str = sprintf('MH:%s',s); % 25 spaces for label name
str(end+1:30) = ' '; % add blanks
for i = 1:length(markerID(:,1))
    idx(i) = strcmp(markerID(i,:),str);
end
indMarker = find(idx == 1,1,'first');
