function [pf, pfu, pSrT, pSrTu] = plotCompfSrT(subj_array,cSubj,indParray,Ulev,pNames,fd,SrTd)
% just plot step freq (f) and stride time (SrT)
% Plot solo, unpulsed, pulsed

figure;
numcols = 2; numrows = 2;
xl = [0.5 length(subj_array)+0.5];

for param = 1:2
    indP = indParray(param);
    indUP = Ulev;
    
    for isubj = 1:length(subj_array) 
        fname = sprintf('HRI%i_metrics.mat',subj_array(isubj));
        data = load(fname); 
        indS = find(data.condArray(:,1) == 0); % solo data
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
                
        % f         
        subplot(numrows,numcols,(param-1)*numcols+1), hold on;
        errorbar(isubj-0.2,nanmean(data.SF(indS)),nanstd(data.SF(indS)),'^','color',cSubj(isubj,:)); % solo
        errorbar(isubj,nanmean(data.SF(indU)),nanstd(data.SF(indU)),'o','color',cSubj(isubj,:)); % unpulsed
        errorbar(isubj+0.2,nanmean(data.SF(indC)),nanstd(data.SF(indC)),'x','color',cSubj(isubj,:)); % pulsed
        ylim([1.3 1.95])
        
        [pfu(isubj,param), stat, test] = compMean(data.SF(indU),fd);
        if pfu(isubj,param) < 0.05
            plot(isubj,nanmean(data.SF(indU))+0.1,'k+'); % denote NR
        end
        
        [pf(isubj,param), stat, test, df] = comp2groups(data.SF(indU),data.SF(indC));
        if pf(isubj,param) < 0.05
            sigstar({[isubj isubj+0.2]});
        end
        
        if isubj == 1 && param == 1
            legend('solo','unpulsed','pulsed');
        elseif isubj == length(subj_array)
%             hline(fd,'k--');
            ylabel('f (hz)'),xlabel('Subj');
        end
        title(pNames{param});
        
        % SrT        
        subplot(numrows,numcols,(param-1)*numcols+2), hold on;
        errorbar(isubj-0.2,nanmean(data.SrT(indS)),nanstd(data.SrT(indS)),'^','color',cSubj(isubj,:)); % solo
        errorbar(isubj-0.2,nanmean(data.SrT(indS)),0.067*nanmean(data.SrT(indS)),'color',cSubj(isubj,:),'linewidth',1.5); % solo BoE
        errorbar(isubj,nanmean(data.SrT(indU)),nanstd(data.SrT(indU)),'o','color',cSubj(isubj,:)); % unpulsed
        errorbar(isubj,nanmean(data.SrT(indU)),0.067*nanmean(data.SrT(indU)),'color',cSubj(isubj,:),'linewidth',1.5); % unpulsed BoE
        errorbar(isubj+0.2,nanmean(data.SrT(indC)),nanstd(data.SrT(indC)),'x','color',cSubj(isubj,:));
        ylim([.9 1.6])
        
        [pSrTu(isubj,param), stat, test] = compMean(data.SrT(indU),SrTd);
        if pSrTu(isubj,param) < 0.05
            plot(isubj,nanmean(data.SrT(indU))+0.1,'k+'); % denote NR
        end    
        
        [pSrT(isubj,param), stat, test, df] = comp2groups(data.SrT(indU),data.SrT(indC));
        if pSrT(isubj,param) < 0.05
            sigstar({[isubj isubj+0.2]});
        end     
                
        if isubj == length(subj_array)
%             hline(SrTd,'k--');
            ylabel('SrT (s)'),xlabel('Subj');
        end        
       
    end
end
a = findobj(gcf,'type','axes');
set(a,'xlim',xl,'box','off','tickdir','out','xtick',1:length(subj_array));