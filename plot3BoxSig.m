function plot3BoxSig(plotind,numrows,numcols,m,xticklab,xl,pval)

subplot(numrows,numcols,plotind)
boxplot(m,'colors','k','symbol','o'); hold on;
box off,set(gca,'xtick',1:length(xticklab),'xticklabel',xticklab,'xlim',xl)
if pval(1) < 0.05
    sigstar({[1 2]});
end
if pval(2) < 0.05
    sigstar({[2 3]});
end
if pval(2) < 0.05
    sigstar({[1 3]});
end