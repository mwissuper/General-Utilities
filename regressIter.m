function [c, rsq, p] = regressIter(y,x)

% y is nx1 dependent variable and x is nx4 matrix of predictors
% Regress iteratively if CI of coeff includes zero
% and remove the predictor from the model. Update output c with zero for
% the removed predictor and update rsq for model.

% Use z to calculate standardized coefficients so can compare their
% relative importance to each other
for i = 1:3
    z(:,i) = zscore(x(:,i));
end
z(:,4) = x(:,4); % ones

% Calculate unstandardized coeff's so they have meaningful physics units
% for mass, damping, stiffness
[b,bint,r,rint,stats] = regress(y,x);
c = b
%% Need to do checks as long as there are predictors that are n.s.

% Check if any predictors are n.s.
for i = 1:length(b)-1 % ignore constant term
    if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
        c(i) = nan; % set coeff to nan so know later which term (i.e. acc, vel, or pos) to drop from x
        c = c
    end
end
% Remove n.s. predictors and regress again if any predictors are n.s.
ind = find(isnan(c)==1,1,'first');
if ~isempty(ind)
    x(:,ind) = 0; % don't delete the column completely so we can keep order of elements of c consistent as m,b,k
end
n = 0;

while ~isempty(ind) && n < length(b) % Keep doing if there are n.s. predictors and 1 or more sig. predictors
    n = n + 1;

    [b,bint,r,rint,stats] = regress(y,x);   
    c = b
    for i = 1:length(b)-1
        if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
            c(i) = nan; 
            c = c
        end
    end
    ind = find(isnan(c)==1,1,'first');
    if ~isempty(ind)
        x(:,ind) = 0;
    end
end

rsq = stats(1);
p = stats(3);

%% Regress one more time to get standardized coefficients
ind = find(isnan(c)==1,1,'first');
if ~isempty(ind)
    z(:,ind) = 0; 
end
[b_st,bint,r,rint,stats] = regress(y,z); 
rsq_st = stats(1);

if rsq_st ~= rsq
    disp('R^2 not equal!')
end