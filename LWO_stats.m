%% Do correlation analysis and also test difference between groups for final lane width data
clear; clc; close all;

% %% Correlations of clin scores with LWO from baseline test before expt -
% % none significant (agrees with case study)
% [metrics, temp] = xlsread('SCI clinical scores.xlsx',1,'C4:S15'); % Skip nan entry for LWO from baseline test
% 
% for i = 2:length(metrics)
%     notnan = find(isnan(metrics(:,i))~=1);
%     [rho(i-1), pval(i-1)] = corr(metrics(notnan,1),metrics(notnan,i));
% end

%% Correlations of clin scores with LWO from Null trial - significance!
clear; clc;
[metrics, temp] = xlsread('SCI clinical scores.xlsx',1,'F3:S14'); % Skip nan entry for LWO from baseline test
[temp, names] = xlsread('SCI clinical scores.xlsx',1,'F2:S2');

for i = 3:length(metrics) % skip col 2 because just one value
    notnan = find(isnan(metrics(:,i))~=1); % skip nan values to calculate correlation
    [Hl(i),Pl(i)] = lillietest(metrics(notnan,i)); % Check if distribution is normal. ALready checked that LWO is normal (col 1 data)
    if Pl(i) < 0.05
        [rho(i), pval(i)] = corr(metrics(notnan,1),metrics(notnan,i),'type','spearman');
    else
        [rho(i), pval(i)] = corr(metrics(notnan,1),metrics(notnan,i));
    end
end

ind = find(pval(3:end) < 0.05) + 2;

for i = 1:length(ind)
    subplot(2,3,i)
    plot(metrics(notnan,1),metrics(notnan,ind(i)),'k.');
    titlename = sprintf('rho = %.2f, p = %.2f',rho(ind(i)),pval(ind(i)));
    title(titlename);
    clear a; box off; axis square; set(gca,'tickdir','out');
    xlabel('LWO width Null trial (m)'),ylabel(names{ind(i)});  
end

% %% Plot for correlation analysis
% % names = temp(1,:);
% [rho,p] = corr(metrics,'rows','complete'); % Ignore nan values
% imagesc(rho);
% colorbar
% set(gca,'xtick',1:length(names),'ytick',1:length(names),'xticklabel',names,'yticklabel',names);
% figure
% corrplot(metrics);

%% Compare 2 groups for Null condition only to see if LWO can distinguish between groups

% Do all stats tests in Matlab to compare results vs. SPSS
% Use lillietest since most comparable to KS test with Lilliefors
% correction in SPSS to determine normality. If normal, use REMANOVA,
% otherwise use Friedman's.

% Within-subject design only! Use R and glmer for mixed models!

[data, temp] = xlsread('LWO final lane widths.xlsx',1,'A2:D73'); % Skip nan entry for LWO from baseline test
group = data(:,1);
cond = data(:,3);
laneWidth = data(:,4);
indisciNull = find(group == 1 & cond == 2);
indniNull = find(group == 2 & cond == 2);

[Hl,P1] = lillietest(laneWidth(indisciNull)); % First test if distrib is normal so know which test to use next
[H2,P2] = lillietest(laneWidth(indniNull)); % First test if distrib is normal so know which test to use next
if P1 < 0.05 | P2 < 0.05
    [P0,H,STATS] = signrank(laneWidth(indisciNull),laneWidth(indniNull));
    W = STATS.signedrank;
else
    [H,P0,CI,STATS] = ttest(laneWidth(indisciNull),laneWidth(indniNull));
    t = STATS.tstat;
end
