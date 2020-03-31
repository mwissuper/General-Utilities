function [c, rsq, p] = regressIter(y,x)

% y is nx1 dependent variable and x is nx3 matrix of predictors
% Regress iteratively if CI of coeff includes zero
% and remove the predictor from the model. Update output c with zero for
% the removed predictor and update rsq for model.

[b,bint,r,rint,stats] = regress(y,x);
c = b;
%% Need to do checks as long as there are predictors that are n.s.

% Check if any predictors are n.s.
for i = 1:length(b)-1 % ignore constant term
    if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
        c(i) = nan; % set coeff to nan so know later which terms (i.e. m, b, or k) were dropped
    end
end
% Remove n.s. predictors and regress again if any predictors are n.s.
ind = find(isnan(c)==1,1,'first');
if ~isempty(ind)
    x(:,ind) = 0;
end
n = 0;

while ~isempty(ind) && n < length(b) % Keep doing if there are n.s. predictors and 1 or more sig. predictors
    n = n + 1;

    [b,bint,r,rint,stats] = regress(y,x);   
    for i = 1:length(b)-1
        if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
            c(i) = nan; 
        end
    end
    ind = find(b==0);
    if ~isempty(ind)
        x(:,ind) = 0;
    end
end

rsq = stats(1);
p = stats(3);