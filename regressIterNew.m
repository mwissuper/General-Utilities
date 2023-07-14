function [c, rsq, p, b_st, rsq_st] = regressIterNew(y,x)

% y is nx1 dependent variable and x is nxm matrix of predictors
% Regress iteratively if CI of coeff includes zero
% and remove the predictor from the model. Update output c with zero for
% the removed predictor and update rsq for model.

% Use z to calculate standardized coefficients so can compare their
% relative importance to each other
s = size(x);
for i = 1:s(2)-1
    z(:,i) = zscore(x(:,i));
end
z(:,i+1) = x(:,end); % ones

% Calculate unstandardized coeff's so they have meaningful physics units
% for mass, damping, stiffness
[b,bint,r,rint,stats] = regress(y,x);
c = b;
blast = b;
%% Need to do checks as long as there are predictors that are n.s.

% Check if any predictors are n.s.
for i = 1:length(b)-1 % ignore constant term
    if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
        c(i) = nan; % set coeff to nan so know later which term (i.e. acc, vel, or pos) to drop from x
%         c = c
        blast(i) = 0;
    end
end
% Remove n.s. predictors and regress again if any predictors are n.s.
ind = find(isnan(c)==1,1,'first');
if ~isempty(ind)
    x(:,ind) = 0; % don't delete the column completely so we can keep order of elements of c consistent as m,b,k
end
n = 0;
flag = 0;
while ~isempty(ind) && length(ind) < s(1) && flag == 0 % Keep doing if there are n.s. predictors and 1 or more sig. predictors
    n = n + 1;

    [b,bint,r,rint,stats] = regress(y,x); 
    if b == blast % same result as last time
        flag = 1;
    end
    blast = b;
    for i = 1:length(b)-1
        if bint(i,1) < 0 && bint(i,2) > 0 || bint(i,1) == 0 || bint(i,2) == 0
            c(i) = nan; 
%             c = c
        end
    end
    ind = find(isnan(c)==1);
    if ~isempty(ind)
        x(:,ind) = 0;
    end
end

rsq = stats(1);
p = stats(3);

%% Regress one more time to get standardized coefficients
ind = find(isnan(c)==1);
if ~isempty(ind)
    z(:,ind) = 0; 
end
[b_st,bint,r,rint,stats] = regress(y,z); 
rsq_st = stats(1);

if abs(rsq_st - rsq) > 1e-5 % sometimes not exactly equal due to digit precision
    disp('R^2 not equal!')
end