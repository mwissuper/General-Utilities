function indMarker = findMarkerIndGen(subj,s,markerID)

% Inputs: string 's' is name of marker and string 'subj' is name of second subject/model in Vicon
% Find which index marker 's' is in MarkerID. Expects length of filename <=
% 30 char's
if isnan(subj) % no subj prefix
    str = s; 
else
    str = sprintf('%s:%s',subj,s); 
end
str(end+1:30) = ' '; % add blanks for file format
for i = 1:length(markerID(:,1))
    idx(i) = strcmp(markerID(i,:),str);
end
indMarker = find(idx == 1,1,'first');
