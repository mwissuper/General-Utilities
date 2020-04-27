function [x, ind] = replaceZeros(x)

% Look for repeated zeros in marker data since Nexus fills gaps with
% zeros. Replace with nans before processing. Find anywhere with value of 0
% repeated at least once. Look at only one component and then replace all
% components with nan's for those samples. Takes in input of nx3 array
% where n is samples of marker data and 3 for x,y,z dir's

% Indices of x that are nan
ind = zeros(size(x(:,2)));

a = find(x(:,2)==0);
b = diff(a);
c = find(b==1); % Indices of x where two zeros in a row
if ~isempty(c)
    for i = 1:3 % all 3 dimensions/components of marker data
        x(a(c),i) = nan;
        x(a(c(end))+1,i) = nan; % Get last zero
    end

    ind(a(c)) = 1;
    ind(a(c(end))+1) = 1; % Get last zero
end