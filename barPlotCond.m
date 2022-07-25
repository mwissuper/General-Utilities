function h = barPlotCond(metric,condLab,titlename,yLabName)

hold on;
n = length(condLab); % number of conditions
h = bar(mean(metric),'linestyle','none');
errorbar(mean(metric),std(metric),'linestyle','none','color','k');
% boxplot(metric,'colors','k','symbol','o')
set(gca,'xticklabel',condLab,'xtick',1:n); box off; set(gca,'tickdir','out');
xlim([0.5 n+0.5]);
title(titlename),ylabel(yLabName); 

