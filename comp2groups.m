function [P0, stat, test] = comp2groups(m1,m2)

% Test normality, then use t-test or non-param test compare 2 paired
% samples
% Two-tailed t-test
[Hl,P1] = lillietest(m1); % First test if distrib is normal so know which test to use next
[H2,P2] = lillietest(m2); % First test if distrib is normal so know which test to use next
if P1 < 0.05 | P2 < 0.05
    if length(m1) == length(m2) 
        test = 'sign rank';
    %     [P0(1),H(1),STATS] = signrank(m1,m2);
    %     stat(1) = STATS.signedrank;
        % Do approximate method to report z val
        [P0(1),H(1),STATS] = signrank(m1,m2,'method','approximate');
    else
        test = 'rank sum';
        [P0(1),H(1),STATS] = ranksum(m1,m2,'method','approximate');
    end
    stat(1) = STATS.zval;
else
    if length(m1) == length(m2) 
        [H,P0,CI,STATS] = ttest(m1,m2);
        test = 't test';
    else
        [H,P0,CI,STATS] = ttest2(m1,m2);
        test = 'Welch t test';
    end
    stat = STATS.tstat;    
end