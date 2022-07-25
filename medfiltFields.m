function S = medfiltFields(S,dim)
%% MEDFILTFIELDS: Median filters all the fields of a structure
%
%   S = medfiltFields(S,dim) applies a 1 dimensional median filter to all
%   the nonscalar fields in S. Each nonscalar field is assumed to be an
%   array, and the filter is applied along dimension DIM of the array.

%   Luke Drnach
%   December 4, 2018

fields = fieldnames(S);
for n = 1:length(fields)
    if ~isscalar(S.(fields{n}))
        S.(fields{n}) = medfilt1(S.(fields{n}),[],[],dim);
    end
end
end