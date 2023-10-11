function ind = selectRefHS(HSss,HSref)

% Select correct HSref events to plot HS events for plotSLST
% inputs are indices of: 
% 1) HS events for foot we want to plot that occurred during 
% SS period, 2) all ref foot HS events3) ref foot HS events during SS period

for i = 1:length(HSss)
    a = find(HSref < HSss(i),1,'last');
    ind(i) = HSref(a);
end

