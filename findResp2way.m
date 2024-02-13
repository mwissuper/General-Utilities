function [resp,pUP,pInt,pP,pph_P] = findResp2way(dv,v1names,v2names,v1ind)
% Find if subj responds to pulses
% Use this code for uninstructed participants where velocity and pulse freq
% vary together. v1names for unpulsed levels. v2names for pulsed levels
% factor 1 names/levels encoded in v1names 
% factor 2 names/levels encoded in v2names 
% dv: dependent var matrix, v1ind: col's of dv 
% that correspond to each level of v1names

%% 2-way REMANOVA
% create table for ranova
c1 = []; c2 = []; c3 = [];
ind = 0;
for j = 1:length(v1names) % num levels of v1
    for k = 1:length(v2names) % num levels of v2 
        c3 = [c3; dv(:,v1ind{j}(k))];
        for i = 1:length(dv(:,1))
            ind = ind + 1;
            c1{ind} = v1names{j};
            c2{ind} = v2names{k};            
        end
    end
end

t = table(c1',c2',c3);

% codes for column names: u = unpulsed, p = pulsed, 1 = level 1, 2 = level
% 2, 3 = level 3
dv = array2table(reshape(t.c3, 5, [])); % dependent variable
dv.Properties.VariableNames = {'u1', 'u2', 'u3', 'p1','p2','p3'};
% create the within-subjects design
withinDesign = table([1 1 1 2 2 2]',[1 2 3 1 2 3]','VariableNames',{'Pulse','Level'});
withinDesign.Pulse = categorical(withinDesign.Pulse);
withinDesign.Level = categorical(withinDesign.Level);
% create repeated measures model
rm = fitrm(dv, 'u1-p3 ~ 1', 'WithinDesign', withinDesign);
% Sphericity...doesn't seem to work: p = 0 for HRI01
sph = mauchly(rm);
% perform anova (remove semicolon to view table generated by ranova)
AT = ranova(rm, 'WithinModel', 'Pulse*Level');

pUP = AT.pValue(2); pV2 = AT.pValue(3); pInt = AT.pValue(4);

%% Follow-up
if pInt < 0.05 % sig interaction of pulsed/unpulsed and velocity level
    resp = 1; 
    pP = nan; pph_P = nan;
    %%  Sig interaction pulse and level, check  at each level whether pulse vs unpulsed are sig diff using EMMEANS
    if pInt < 0 
        % Need estimated marginal means to follow up which interaction
        % effects are sig.
    else
        pphInt = nan*ones(3,1);
    end
else % no sig interaction, check main effect of pulsed/unpulsed
    pphInt = nan;
    if pUP < 0.05 
        resp = 1;
        % It will be obvious which has a higher mean since 2 levels of this
        % factor
    else
        resp = 0;
        pP = nan; pph_P = nan;
    end
end

