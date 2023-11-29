function [resp,pU,pP] = findRespUnbal(dv,v1names,v1ind)
% Multi-step stats tests to determine if subj responds to pulses
% factor 1 names/levels encoded in v1names and 
% dv: dependent var matrix, v1ind: col's of dv corresponding to levels of
% factor 1
% Assumes different num levels for U as P condition, so can't do 2-way
% REMANOVA

%% 1-way REMANOVA to determine if w is constant in unpulsed. Assume 3 levels for unpulsed.
[sph_p,ranovatbl,pRA,F,df1,df2,posthoctbl,pph] = run_ranova_3lev(dv(:,v1ind{1}));
pU = pRA;
if pU >= 0.05
    % 1-way REMANOVA to dtermine if w is constant in pulsed.
    [sph_p,ranovatbl,pRA,F,df1,df2,posthoctbl,pph] = run_ranova_3lev(dv(:,v1ind{2}));
    pP = pRA;
    if pP < 0.05
        resp = 1;
    else
        resp = 0;
    end
else
    pP = nan;
    resp = nan;
end

