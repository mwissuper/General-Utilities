function h = plotGaitError(subj_array,cSubj,indParray,Ulev,vd,meanV,pv,pvu,fd,meanSF,pf,pfu,Ld,meanSL,pL,pLu,wd,meanW,pw,pwu)
% Plot error in gait metrics for each participant in pulsed condition
% relative to unpulsed condition

numcols = 4; numrows = 1;
xl = [0.5 length(subj_array)+0.5];
pNames = {'dP','va'};
h = figure;
offset = [-0.2 0.2];
for param = 1:2
    indC = indParray(param);
    for isubj = 1:length(subj_array) 

        %% v    
        yl = [-1 1]*.017;
        dv = abs(meanV(isubj,indC)-vd) - abs(meanV(isubj,Ulev)-vd);
        subplot(numrows,numcols,1), hold on;
        bar(isubj+offset(param),dv,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pv(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pvu(isubj,param) < 0.05
                plot(isubj+offset(param),dv+sign(dv)*yl(2)*0.1,'k+');
            end
        end
        ylabel('dvError (m/s)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',yl)
        
        %% f   
        yl = [-1 1]*.16;
        df = abs(meanSF(isubj,indC)-fd) - abs(meanSF(isubj,Ulev)-fd);
        subplot(numrows,numcols,2), hold on;
        bar(isubj+offset(param),df,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pf(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pfu(isubj,param) < 0.05
                plot(isubj+offset(param),df+sign(df)*yl(2)*0.1,'k+');
            end
        end
        ylabel('dfError (hz)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',yl)
        
        %% L         
        yl = [-1 1]*.05;
        dL = abs(meanSL(isubj,indC)-Ld) - abs(meanSL(isubj,Ulev)-Ld);
        subplot(numrows,numcols,3), hold on;
        bar(isubj+offset(param),dL,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pL(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pLu(isubj,param) < 0.05
                plot(isubj+offset(param),dL+sign(dL)*yl(2)*0.1,'k+');
            end
        end
        ylabel('dLError (m)'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',yl)
        
        %% w         
        yl = [-1 1]*.3;
        dw = abs(meanW(isubj,indC)-wd) - abs(meanW(isubj,Ulev)-wd);
        subplot(numrows,numcols,4), hold on;
        bar(isubj+offset(param),df,'facecolor',cSubj(isubj,:),'edgecolor','k','barwidth',0.4);
        if pw(isubj,param) < 0.05
            sigstar({isubj+offset(param)});
        else
            if pwu(isubj,param) < 0.05
                plot(isubj+offset(param),dw+sign(dw)*yl(2)*0.1,'k+');
            end
        end
        ylabel('dwError'),xlabel('Subj'),set(gca,'xtick',1:length(subj_array),'xlim',xl,'ylim',yl)
   
    end      
end