function [P0, stat, test] = comp2groups(m1,m2)

% Test normality, then use t-test or sign-rank test compare 2 paired
% samples
% Two-tailed t-test
[Hl,P1] = lillietest(m1); % First test if distrib is normal so know which test to use next
[H2,P2] = lillietest(m2); % First test if distrib is normal so know which test to use next
if P1 < 0.05 | P2 < 0.05
    test = 'sign rank';
%     [P0(1),H(1),STATS] = signrank(m1,m2);
%     stat(1) = STATS.signedrank;
    % Do approximate method to report z val
    [P0(1),H(1),STATS] = signrank(m1,m2,'method','approximate');
    stat(1) = STATS.zval;
else
    [H,P0,CI,STATS] = ttest(m1,m2);
    test = 't test';
    stat = STATS.tstat;
end