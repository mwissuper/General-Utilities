function p = repDesignStats(data,num_conds,refcol,condNames)

% Function to do stats analysis for all metrics for 5 cond's

% Do Bonferroni correction of number of pairwise comparisons vs. Null
if num_conds == 3
    bonf = 3;
elseif num_conds == 5
    bonf = 4;
end
% Remove nan's if missing subj
[i,j] = find(isnan(data)==1);
data(i,:) = [];

compcols = 1:num_conds;
compcols(refcol) = []; % Remove the reference condition

% KS test
for i = 1:num_conds
    [hKS(i), pKS(i)] = kstest(data(:,i));
end

if isempty(find(hKS==1,1,'first')) == 0 % Normality violated for some condition
    % Check omnibus effect with friedman's test
    pFriedman = friedman(data,1) % Print it out
    if pFriedman < 0.05
        % Do Wilcoxon's test for each condition vs. baseline, adjust p to
        % account for multiple comparisons
        n = 0;
        for i = [1 2 4 5]
            n = n + 1;
            [p(n),h,stats] = signrank(data(:,i),data(:,3),'method','approximate'); % Bonferroni correction
            if p(n)*bonf < 0.05 % Bonferroni correction for 6 comparisons
                if median(data(:,3)) < median(data(:,i))
                    sprintf('Null < %s p = %.3f, z = %.2f',condNames{i},p(n)*bonf,stats.zval)
                else
                    sprintf('Null > %s p = %.3f, z = %.2f',condNames{i},p(n)*bonf,stats.zval)
                end
            end
        end
    end
else % Do REMANOVA
    % Check sphericity
end

close all;