function [sph_p,ranovatbl,pRA,F,df1,df2,posthoctbl,pph] = run_ranova_3lev(m)

% m is metric of interest, assume 3 levels of conditions. m is n x 3 matrix
% where n is number of subjects 

t = table(m(:,1),m(:,2),m(:,3),'VariableNames',{'cond1','cond2','cond3'});
withinDesign = table([1 2 3]','VariableNames',{'Condition'});
withinDesign.Condition = categorical(withinDesign.Condition);
rm = fitrm(t,'cond1-cond3 ~ 1','WithinDesign',withinDesign);
sph = mauchly(rm);
sph_p = sph.pValue;
ranovatbl = ranova(rm);

% Get F and df
F = ranovatbl.F(1);
df1 = ranovatbl.DF(1);
df2 = ranovatbl.DF(2);

% Use appropriate p value depending on sphericity
if sph.pValue > 0.05 
    pRA = ranovatbl.pValue(1);
else % use gg p, calculate df's
    pRA = ranovatbl.pValueGG(1);
    etbl = epsilon(rm);
    e = etbl.GreenhouseGeisser;
    df1 = df1*e;
    df2 = df2*e;
end

% Do follow-up if p is sig
if pRA < 0.05
    posthoctbl = multcompare(rm,'Condition','ComparisonType','bonferroni');
    pph = posthoctbl.pValue([1 2 4]); % Extract out p12, p13, p23
else 
    posthoctbl = nan;
    pph = nan*ones(1,3);
end