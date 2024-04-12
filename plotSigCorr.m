function [r, p] = plotSigCorr(x,y,numrows,numcols,plotind)

[p, r, test] = corr2vars(x,y);
% if p < 0.05
    subplot(numrows,numcols,plotind)
    c = polyfit(x,y,1);
    plot(x,polyval(c,x),'k-'); 
% end
tname = sprintf('p: %.2f, r: %.3f',p,r); title(tname);