function b = plotBoxSig(plotind,numrows,numcols,m,xticklab,xl,pval)

subplot(numrows,numcols,plotind)
b = boxplot(m,'colors','k','symbol','o'); hold on;
box off,set(gca,'xtick',1:length(xticklab),'xticklabel',xticklab,'xlim',xl)
if length(pval) == 3
    pmatrix = [1 2;1 3;2 3];
elseif length(pval) == 4 
%     pmatrix = [1 2;2 3;2 4];
end
if ~isnan(pval)
    for i = 1:length(pval)
        if pval(i) < 0.05
            sigstar({pmatrix(i,:)});
        end
    end
end