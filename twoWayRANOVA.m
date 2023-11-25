function twoWayRANOVA(dv,v1names,v2names,v1ind)
% 2-way rep measures ANOVA with factor 1 levels encoded in v1 names and 
% factor 2 levels encoded in v2 names. 
% dv: dependent var matrix, v1ind: col's of dv
% that correspond to each level of v1
% length of v1ind encodes levels of factor 2

% create col's of v1 and v2
for j = 1:length(v1ind) % num levels of v2
    ind = ((j-1)*length(dv,1)+1):(j*length(dv,1));
    c1(ind) = 1:length(dv,1); % subj id
    for i = 1:length(dv,1)
        c2{i} = v1names;
    end
end