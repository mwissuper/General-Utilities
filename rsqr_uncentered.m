function ursqr = rsqr_uncentered(data,data_rec)
% This function calculates the uncetered correlation coefficient using "Cluster" method.  
%
% Syntax:   r_sqr = rsqr_uncentered(data,data_rec)
%
% Input:
% data      Array   matrix of observed data  (e.g., data = [observations x channels])
% data_rec  Array   matrix of reconstructed/predicted data (e.g., data_rec = [observations x channels])
%
% Output:
% ursqr     Array   row vector of uncentered correlation coefficients
%       
% Created: May 24, 2006 (Gelsy Torres-Oviedo)
% Last Modified: July 10, 2006 (Torrence Welch)
% Last Modification: fix ursqr calculation
% Last Modified: June 7, 2016 (J. Lucas McKay)
% Last Modification:
%   remove transpose to expect data in columns
%   remove empty values

% Zar book method p. 334
nvars = size(data,2);
ursqr = nan(1,nvars);
for i = 1:nvars
    X = [data(:,i) data_rec(:,i)];
    X(isnan(X(:,1))|isnan(X(:,2)),:) = [];    
    ursqr(i) = sum(prod(X,2))^2 / (sum(X(:,1).^2)*sum(X(:,2).^2)); %regression sum of squares/total sum of squares
end

end