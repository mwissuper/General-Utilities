function [P0, stat, test] = compMean(a,m)

% Test normality, then use t-test or sign-rank test compare sample a
% against mean m
% Two-tailed t-test
[Hl,P1] = lillietest(a); % First test if distrib is normal so know which test to use next

if P1 < 0.05
    [P0(1),H(1),STATS] = signrank(a,m);
    test = 'sign rank';
    stat(1) = STATS.signedrank;
    % Also do approximate method to report z val
    [P0(2),H(2),STATS] = signrank(a,m,'method','approximate');
    stat(2) = STATS.zval;
else
    [H,P0,CI,STATS] = ttest(a,m);
    test = 't test';
    stat = STATS.tstat;
end