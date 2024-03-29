function [pV1,pV2,pInt,pphInt] = twoWayRANOVA(v,v1names,v2names,vInd)
% 2-way rep measures ANOVA with factor 1 names/levels encoded in v1names and 
% factor 2 names/levels encoded in v2names. 
% dv: dependent var matrix, v1ind: col's of dv 
% that correspond to each level of v1

c1 = []; c2 = []; c3 = []; s = []; f1 = []; f2 = [];
ind = 0;
% create table for ranova
for j = 1:length(v1names) % num levels of v1
    for k = 1:length(v2names) % num levels of v2 
        s = [s; (1:length(v(:,1)))'];
        f1 = [f1; j*ones(length(v(:,1)),1)];
        f2 = [f2; k*ones(length(v(:,1)),1)];
        c3 = [c3; v(:,vInd{j}(k))];
        for i = 1:length(v(:,1))
            ind = ind + 1;
            c1{ind} = v1names{j};
            c2{ind} = v2names{k};            
        end
    end
end

t = table(s,f1,f2,c3);

%% Use function from Matlab Central
% t.Properties.VariableNames = {'cond','level','w'};
% stats = rm_anova2(c3,s,f1,f2,{'cond','level'})
% pV1 = stats{2,5};
% pV2 = stats{3,5};
% pInt = stats{4,5};

%% Adapted this code://www.mathworks.com/matlabcentral/answers/886394-repeated-measures-anova-with-two-variables?s_tid=mlc_ans_email_ques

% codes for column names: u = unpulsed, p = pulsed, 1 = level 1, 2 = level
% 2, 3 = level 3
% organize the data in a matrix (1 row per participant)
dv = array2table(reshape(t.c3, length(v(:,1)), [])); % dependent variable is 'Value'
dv.Properties.VariableNames = {'u1','u2','u3','p1','p2','p3'};
% create the within-subjects design
withinDesign = table([1 1 1 2 2 2]',[1 2 3 1 2 3]','VariableNames',{'cond','level'});
withinDesign.cond = categorical(withinDesign.cond);
withinDesign.level = categorical(withinDesign.level);
% create repeated measures model
rm = fitrm(dv, 'u1-p3 ~ 1', 'WithinDesign', withinDesign);
% Sphericity...doesn't seem to work: p = 0 for HRI01
sph = mauchly(rm);
% perform anova (remove semicolon to view table generated by ranova)
AT = ranova(rm, 'WithinModel', 'cond*level');

pV1 = AT.pValue(3); pV2 = AT.pValue(5); pInt = AT.pValue(7);

% % post hoc comparisons between all columns (2 factors merged) --> different
% % result than EMMEANS in SPSS, don't trust this code.
% withinDesign2 = table([1 2 3 4 5 6]', 'VariableNames', {'IV2'});
% withinDesign2.IV2 = categorical(withinDesign2.IV2);
% rm2 = fitrm(dv, 'u1-p3 ~ 1', 'WithinDesign', withinDesign2);
% multcompare(rm2, 'IV2')

% %%  Sig interaction pulse and level, check  at each level whether pulse vs unpulsed are sig diff
% if pInt < 0.05 

% else 
    pphInt = nan*ones(3,1);
% end

