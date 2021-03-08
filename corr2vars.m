function [p, r, test] = corr2vars(v1,v2)

% Test normality, then use t-test or sign-rank test compare 2 paired
% samples
% Two-tailed t-test
[Hl,P1] = lillietest(v1); % First test if distrib is normal so know which test to use next
[H2,P2] = lillietest(v2); % First test if distrib is normal so know which test to use next
if P1 < 0.05 | P2 < 0.05
    test = 'Spearman';
    [r,p] = corr(v1,v2,'type','spearman');
else % Pearson's corr
    test = 'Pearson';
    [r,p] = corr(v1,v2);
end