function [h, impdf, worsedf] = plotGaitError(subj_array,cSubj,indParray,Ulev,vd,meanV,pv,pvu,fd,meanSF,pf,pfu,Ld,meanSL,pL,pLu,wd,meanW,pw,pwu)
% Plot error in gait metrics for each participant in pulsed condition
% relative to unpulsed condition. Need yl if want to label NR with custom
% marker a set distance above the bar

numcols = 4; numrows = 1;
xl = [0.5 length(subj_array)+0.5];
pNames = {'dP','va'};
h = figure;
offset = [-0.2 0.2]; % x offset for bars
ylv = [-1 1]*.045;
ylf = [-1 1]*.17;
ylL = [-1 1]*.08;
ylw = [-1 1]*.15;


        
for param = 1:2
    indC = indParray(param);
    
    % Print how many improved error per pulse shape
    dfAll = abs(meanSF(:,indC)-fd) - abs(meanSF(:,Ulev)-fd);
    impdf(param) = length(find(dfAll < 0 & pf(:,param) < 0.05)); 
    worsedf(param) = length(find(dfAll > 0 & pf(:,param) < 0.05)); 
%     numImpdf(param) = length(find(dfAll < 0))
%     
    % Not sure why this code for yl didn't work. Anyway, want uniform ylim
    % across all pulsed cond's. Maybe didn't work because adjusts yl per
    % pulse shape!
%     ylv = [-1.1 1.1]*max(max(abs( abs(meanV(:,indC)-vd) - abs(meanV(:,Ulev)-vd) )));
%     ylf = [-1.1 1.1]*max(max(abs( abs(meanSF(:,indC)-fd) - abs(meanSF(:,Ulev)-fd) )));
%     ylL = [-1.1 1.1]*max(max(abs( abs(meanSL(:,indC)-Ld) - abs(meanSL(:,Ulev)-Ld) )));
%     ylw = [-1.1 1.1]*max(max(abs( abs(meanW(:,indC)-wd) - abs(meanW(:,Ulev)-wd) )));
    
    for isubj = 1:length(subj_array) 

        %% v            
        dv = abs(meanV(isubj,indC)-vd) - abs(meanV(isubj,Ulev)-vd);
        subplot(numrows,numcols,1), hold on;
        bar(isubj+offset(param),dv,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pv(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pvu(isubj,param) < 0.05
                plot(isubj+offset(param),dv+sign(dv)*ylv(2)*0.05,'k+');
            end
        end
        ylabel('dvError (m/s)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',ylv)
        
        %% f   
        df = abs(meanSF(isubj,indC)-fd) - abs(meanSF(isubj,Ulev)-fd);
        subplot(numrows,numcols,2), hold on;
        bar(isubj+offset(param),df,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pf(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pfu(isubj,param) < 0.05
                plot(isubj+offset(param),df+sign(df)*ylf(2)*0.05,'k+');
            end
        end
        ylabel('dfError (hz)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',ylf)        
        
        %% L         
        dL = abs(meanSL(isubj,indC)-Ld) - abs(meanSL(isubj,Ulev)-Ld);
        subplot(numrows,numcols,3), hold on;
        bar(isubj+offset(param),dL,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pL(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pLu(isubj,param) < 0.05
                plot(isubj+offset(param),dL+sign(dL)*ylL(2)*0.05,'k+');
            end
        end
        ylabel('dLError (m)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',ylL)
        
        %% w         
        dw = abs(meanW(isubj,indC)-wd) - abs(meanW(isubj,Ulev)-wd);
        subplot(numrows,numcols,4), hold on;
        bar(isubj+offset(param),df,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pw(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pwu(isubj,param) < 0.05
                plot(isubj+offset(param),dw+sign(dw)*ylw(2)*0.05,'k+');
            end
        end
        ylabel('dwError'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',ylw)
   
    end      
end
impdf = impdf'; worsedf = worsedf';