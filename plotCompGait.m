function [pv, pvu, pf, pfu, pL, pLu, pw, pwu] = plotCompGait(subj_array,cSubj,indParray,Ulev,pNames,vd,fd,Ld,wd)

figure;
numcols = 4; numrows = 1;
xl = [0.5 length(subj_array)+0.5];

for isubj = 1:length(subj_array) % Load data once per subj
    fname = sprintf('HRI%i_metrics.mat',subj_array(isubj));
    data = load(fname); 
    w = data.SF./data.SL;
    indU = find(data.condArray(:,1) == Ulev);
    % Remove first trial if it's an outlier and there are enough valid
    % trials without it
    if length(indU) == 5 && isempty(find(isnan(data.SF(indU))==1,1,'first'))
        tf = isoutlier(data.SF(indU));
        if tf(1) == 1
            indU(1) = [];
        end
    end
    
    for param = 1:length(indParray)
        indP = find(data.condArray(:,1) == indParray(param));
        % Remove first trial if it's an outlier and there are enough valid
        % trials without it
        if length(indP) == 5 && isempty(find(isnan(data.SF(indP))==1,1,'first'))
            tf = isoutlier(data.SF(indP));
            if tf(1) == 1
                indP(1) = [];
            end
        end
        
        % v         
        subplot(numrows,numcols,1), hold on;
        if param == 1
            errorbar(isubj-0.2,nanmean(data.speed(indU)),nanstd(data.speed(indU)),'o','color',cSubj(isubj,:)); % unpulsed
            hold on;
            errorbar(isubj,nanmean(data.speed(indP)),nanstd(data.speed(indP)),'x','color',cSubj(isubj,:)); % dP
        else
            errorbar(isubj+0.2,nanmean(data.speed(indP)),nanstd(data.speed(indP)),'^','color',cSubj(isubj,:)); %va
        end        
        [pv(isubj,param), stat, test, df] = comp2groups(data.speed(indU),data.speed(indP));
        if pv(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2*(param-1)]});
        end
        [pvu(isubj,param), stat, test] = compMean(data.speed(indU),vd); % If sig, then NR
        if isubj == length(subj_array) && param == 2
            %hline(vd,'k--');
            ylabel('v (m/s)'),xlabel('Subj'),legend('Unpulsed','dp','va');
        end
        
        % f         
        subplot(numrows,numcols,2), hold on;
        if param == 1
            errorbar(isubj-0.2,nanmean(data.SF(indU)),nanstd(data.SF(indU)),'o','color',cSubj(isubj,:));
            errorbar(isubj,nanmean(data.SF(indP)),nanstd(data.SF(indP)),'x','color',cSubj(isubj,:));
        else
            errorbar(isubj+0.2,nanmean(data.SF(indP)),nanstd(data.SF(indP)),'^','color',cSubj(isubj,:));
        end
        [pf(isubj,param), stat, test, df] = comp2groups(data.SF(indU),data.SF(indP));
        if pf(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2*(param-1)]});
        end
        [pfu(isubj,param), stat, test] = compMean(data.SF(indU),fd); 
        if isubj == length(subj_array) && param == 2
            %hline(fd,'k--');
            ylabel('f (hz)'),xlabel('Subj');
        end
        
%         % SrT        
%         subplot(numrows,numcols,(param-1)*numcols+3), hold on;
%         errorbar(isubj-0.2,nanmean(data.SrT(indU)),nanstd(data.SrT(indU)),'o','color',cSubj(isubj,:));
%         errorbar(isubj-0.2,nanmean(data.SrT(indU)),0.067*nanmean(data.SrT(indU)),'o','color',cSubj(isubj,:),'linewidth',1.5); % BoE
%         errorbar(isubj+0.2,nanmean(data.SrT(indP)),nanstd(data.SrT(indP)),'x','color',cSubj(isubj,:));
%         [pSrT(isubj,param), stat, test, df] = comp2groups(data.SrT(indU),data.SrT(indP));
%         if pSrT(isubj,param) < 0.05
%             sigstar({[isubj-0.2 isubj+0.2]});
%         end
%         [pSrTu(isubj,param), stat, test] = compMean(data.SrT(indU),SrTd); % plot later, just test here
%         if isubj == length(subj_array)
%             hline(SrTd,'k--');
%             ylabel('SrT (s)'),xlabel('Subj');
%         end
        
        % L        
%         subplot(numrows,numcols,(param-1)*numcols+3), hold on;
        subplot(numrows,numcols,3), hold on;
        if param == 1
            errorbar(isubj-0.2,nanmean(data.SL(indU)),nanstd(data.SL(indU)),'o','color',cSubj(isubj,:));
            errorbar(isubj,nanmean(data.SL(indP)),nanstd(data.SL(indP)),'x','color',cSubj(isubj,:));
        else
            errorbar(isubj+0.2,nanmean(data.SL(indP)),nanstd(data.SL(indP)),'^','color',cSubj(isubj,:));
        end
        [pL(isubj,param), stat, test, df] = comp2groups(data.SL(indU),data.SL(indP));
        if pL(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2*(param-1)]});
        end
        [pLu(isubj,param), stat, test] = compMean(data.SL(indU),Ld); % plot later, just test here
        if isubj == length(subj_array) && param == 2
            %hline(Ld,'k--');
            ylabel('L (m)'),xlabel('Subj');
        end
        
        % w         
        subplot(numrows,numcols,4), hold on;
        if param == 1
            errorbar(isubj-0.2,nanmean(w(indU)),nanstd(w(indU)),'o','color',cSubj(isubj,:));
            errorbar(isubj,nanmean(w(indP)),nanstd(w(indP)),'x','color',cSubj(isubj,:));
        else
            errorbar(isubj+0.2,nanmean(w(indP)),nanstd(w(indP)),'^','color',cSubj(isubj,:));
        end
        [pw(isubj,param), stat, test, df] = comp2groups(w(indU),w(indP));
        if pw(isubj,param) < 0.05
            sigstar({[isubj-0.2 isubj+0.2*(param-1)]});
        end
        [pwu(isubj,param), stat, test] = compMean(w(indU),wd); % plot later, just test here
        if isubj == length(subj_array) && param == 2
            %hline(wd,'k--');
            ylabel('w'),xlabel('Subj');
        end
    end
end

a = findobj(gcf,'type','axes');
set(a,'xlim',xl,'box','off','tickdir','out','xtick',1:length(subj_array));