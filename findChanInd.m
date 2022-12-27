function indChan = findChanInd(s,chanID)

% Take input string s of name of channel and find which index it is in
% chanID. Currently using for EMG data.

str = s; % 24 spaces for label name
str(end+1:24) = ' ';
for i = 1:length(chanID(:,1))
    idx(i) = strcmp(chanID(i,:),str);
end
indChan = find(idx == 1,1,'first');