function [resp,pU,pP] = findRespUnbal(DV,IVnd)
% Multi-step stats tests to determine if subj responds to pulses
% Unbalanced design so can't compare unpulsed vs pulsed directly
% dv: dependent var matrix; row = trial, col = level of factor
% v1ind: struct where each cell of v1ind contains indices of columns of dv matrix that correspond to each indep var.
% Assumes different num levels for U as P condition, so can't do 2-way
% REMANOVA

[sph_p,ranovatbl,pRA,F,df1,df2,posthoctbl,pph] = run_ranova_3lev(DV(:,IVnd{1})); % 1-way REMANOVA to determine if w is constant in unpulsed (IV 1). Assume 3 levels for unpulsed.
pU = pRA;
[sph_p,ranovatbl,pRA,F,df1,df2,posthoctbl,pph] = run_ranova_3lev(DV(:,IVnd{2})); % 1-way REMANOVA to determine if w is constant in pulsed (IV 2).
pP = pRA;
if pU >= 0.05 % w is constant across speeds in unpulsed    
    if pP < 0.05
        resp = 1;
    else
        resp = 0;
    end
else 
    if pP < 0.05
        resp = nan; % both unpulsed and pulsed show changes in w, can'
    else
        resp = 0;
    end
end

