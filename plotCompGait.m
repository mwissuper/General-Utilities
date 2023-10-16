function [pv, pvu, pf, pfu, pL, pLu, pw, pwu] = plotCompGait(subj_array,cSubj,indParray,Ulev,pNames,vd,fd,Ld,wd)

figure;
numcols = 4; numrows = 2;
xl = [0.5 length(subj_array)+0.5];

for param = 1:2
    indP = indParray(param);
    indUP = Ulev;
    
    for isubj = 1:length(subj_array) 
        fname = sprintf('HRI%i_metrics.mat',subj_array(isubj));
        data = load(fname); 
        w = data.SF./data.SL; wp = data.fp/data.Lp;
        indU = find(data.condArray(:,1) == indUP);
        if length(indU) == 5 && isempty(find(isnan(data.SF(indU))==1,1,'first'))
            tf = isoutlier(data.SF(indU));
            if tf(1) == 1
                indU(1) = [];
            end
        end
        indC = find(data.condArray(:,1) == indP);
        if length(indC) == 5 && isempty(find(isnan(data.SF(indC))==1,1,'first'))
            tf = isoutlier(data.SF(indC));
            if tf(1) == 1
                indC(1) = [];
            end
        end
        
        % v         
        subplot(numrows,numcols,(param-1)*numcols+1), hold on;
        errorbar(isubj-0.2,nanmean(data.speed(indU)),nanstd(data.speed(indU)),'o','color',cSubj(isubj,:));
        errorbar(isubj+0.2,nanmean(data.speed(indC)),nanstd(data.speed(indC)),'x','color',cSubj(isubj,:));
        [pv(isubj,param), stat, test, df] = comp2groups(data.speed(indU),data.speed(indC));
        if pv(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2]});
        end
        [pvu(isubj,param), stat, test] = compMean(data.speed(indU),data.vp); % plot later, just test here
        if isubj == length(subj_array)
            hline(vd,'k--');
            ylabel('v (m/s)'),xlabel('Subj'),legend('Unpulsed','Pulsed');
        end
        title(pNames{param});
        
        % f         
        subplot(numrows,numcols,(param-1)*numcols+2), hold on;
        errorbar(isubj-0.2,nanmean(data.SF(indU)),nanstd(data.SF(indU)),'o','color',cSubj(isubj,:));
        errorbar(isubj+0.2,nanmean(data.SF(indC)),nanstd(data.SF(indC)),'x','color',cSubj(isubj,:));
        [pf(isubj,param), stat, test, df] = comp2groups(data.SF(indU),data.SF(indC));
        if pf(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2]});
        end
        [pfu(isubj,param), stat, test] = compMean(data.SF(indU),data.fp); % plot later, just test here
        if isubj == length(subj_array)
            hline(fd,'k--');
            ylabel('f (hz)'),xlabel('Subj');
        end
        
        % L        
        subplot(numrows,numcols,(param-1)*numcols+3), hold on;
        errorbar(isubj-0.2,nanmean(data.SL(indU)),nanstd(data.SL(indU)),'o','color',cSubj(isubj,:));
        errorbar(isubj+0.2,nanmean(data.SL(indC)),nanstd(data.SL(indC)),'x','color',cSubj(isubj,:));
        [pL(isubj,param), stat, test, df] = comp2groups(data.SL(indU),data.SL(indC));
        if pL(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2]});
        end
        [pLu(isubj,param), stat, test] = compMean(data.SL(indU),data.Lp); % plot later, just test here
        if isubj == length(subj_array)
            hline(Ld,'k--');
            ylabel('L (m)'),xlabel('Subj');
        end
        
        % w         
        subplot(numrows,numcols,(param-1)*numcols+4), hold on;
        errorbar(isubj-0.2,nanmean(w(indU)),nanstd(w(indU)),'o','color',cSubj(isubj,:));
        errorbar(isubj+0.2,nanmean(w(indC)),nanstd(w(indC)),'x','color',cSubj(isubj,:));
        [pw(isubj,param), stat, test, df] = comp2groups(w(indU),w(indC));
        if pw(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2]});
        end
        [pwu(isubj,param), stat, test] = compMean(w(indU),wp); % plot later, just test here
        if isubj == length(subj_array)
            hline(wd,'k--');
            ylabel('w'),xlabel('Subj');
        end
    end
end
a = findobj(gcf,'type','axes');
set(a,'xlim',xl,'box','off','tickdir','out','xtick',1:length(subj_array));