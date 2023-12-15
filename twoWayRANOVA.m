function [pV1,pV2,pInt,pphUP,pphInt] = twoWayRANOVA(dv,v1names,v2names,v1ind)
% 2-way rep measures ANOVA with factor 1 names/levels encoded in v1names and 
% factor 2 names/levels encoded in v2names. 
% dv: dependent var matrix, v1ind: col's of dv 
% that correspond to each level of v1

c1 = []; c2 = []; c3 = [];
ind = 0;
% create table for ranova
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

pV1 = AT.pValue(2); pV2 = AT.pValue(3); pInt = AT.pValue(4);

%% Sig diff between levels. Follow up with 1-way REMANOVA of level for unpulsed condition
% --> need to work on this code
if pV2 < 0.05
    posthoctbl = multcompare(rm,'Level','ComparisonType','bonferroni'); % This compares across 
    pphV2 = posthoctbl.pValue([1 2 4]); % Extract out p12, p13, p23
else
    pphV2 = nan*ones(3,1);
end

%%  Sig interaction pulse and level, check  at each level whether pulse vs unpulsed are sig diff
if pInt < 0 
%     posthoctbl =
%     multcompare(rm,'Pulse:Level','ComparisonType','bonferroni'); % can't
%     use this code
    pphInt = posthoctbl.pValue([1 2 4]); % Extract out p12, p13, p23
else
    pphInt = nan*ones(3,1);
end
