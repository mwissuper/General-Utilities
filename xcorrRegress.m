function [c, rsq] = xcorrRegress(y,x)

% y is nx1 dependent variable and x is nx3 matrix of predictors
% First find lag that maximizes xcorr  with y for each predictor, then use
% that lag for regression. Regress iteratively if CI of coeff includes zero
% and remove the predictor from the model. Update output c with zero for
% the removed predictor and update rsq for model.

for i = 1:length(x(1,:))
    d = finddelay(y,x(:,i));
    if d < 0
        xShift = x((1-d):end,i);
        y2 = y(1:end+d);
    elseif d > 0
        xShift = x(1:end-d);
        y2 = y((1+d):end);
    else
        xShift = x(:,i);
        y2 = y;
    end
end




rsq = stats(1);
p = stats(3);