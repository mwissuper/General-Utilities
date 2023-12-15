function [P0, stat, test, df] = comp2groups(m1,m2,dir)

% Test normality, then use t-test or non-param test compare 2 paired
% samples with 2-tailed t-test
if dir == 0
    tail = 'left'; % m1 < m2
elseif dir == 1
    tail = 'right';
else
    tail = 'both';
end

[Hl,P1] = lillietest(m1); % First test if distrib is normal so know which test to use next
[H2,P2] = lillietest(m2); % First test if distrib is normal so know which test to use next
if P1 < 0.05 | P2 < 0.05
    if length(m1) == length(m2) 
        test = 'sign rank';
    %     [P0(1),H(1),STATS] = signrank(m1,m2);
    %     stat(1) = STATS.signedrank;
        % Do approximate method to report z val
        [P0(1),H(1),STATS] = signrank(m1,m2,'method','approximate','tail',tail); % Wilcoxon sign-rank test, T statistic
    else
        test = 'rank sum';
        [P0(1),H(1),STATS] = ranksum(m1,m2,'method','approximate','tail',tail); % Mann-Whitney U test, U statistic
    end
    stat(1) = STATS.zval; % no df required
    df = nan;
else
    if length(m1) == length(m2) 
        [H,P0,CI,STATS] = ttest(m1,m2,'tail',tail); % Student's t test, t statistic
        test = 't test';
    else
        [H,P0,CI,STATS] = ttest2(m1,m2,'tail',tail); % Welch's t test, t statistic
        test = 'Welch t test';
    end
    stat = STATS.tstat; 
    df = STATS.df;
end