function rsqr = rsqr(data,data_rec)
% calculates coefficient of determination (R-squared) according to
% methodology suggested by the Mathworks
%
% https://www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
% accessed 2017 07 07
%
% Created: June 7, 2017 (J. Lucas McKay)
% last edit: 3/20/19: corrected an error in indexing the STATS structure
% that failed if the number of data columns was >1.

[nrows,ncols] = size(data);
rsqr = nan(1,ncols);

% method = categorical(cellstr('fitlm'));
method = categorical(cellstr('regress'));

if method == 'fitlm'
    for i = 1:ncols
        X = data_rec(:,i);
        y = data(:,i);
        mdl = fitlm(X,y);
        rsqr(i) = mdl.Rsquared.Ordinary;
    end
elseif method == 'regress'
    for i = 1:ncols
        X = [data_rec(:,i) ones(size(data_rec(:,i)))];
        y = data(:,i);
        [B,~,~,~,stats] = regress(y,X);
        rsqr(i) = stats(1);
    end
end

% the following instances demonstrate that the regress method is
% approximatley 15x faster than the fitlm method.
if false
    f = [];
    r = [];
    for j = 1:1000
        tic
        X = data_rec(:,i);
        y = data(:,i);
        mdl = fitlm(X,y);
        rsqr(i) = mdl.Rsquared.Ordinary;
        f(end+1) = toc;
        
        tic
        X = [data_rec(:,i) ones(size(data_rec(:,i)))];
        y = data(:,i);
        [B,~,~,~,stats] = regress(y,X);
        rsqr(i) = stats(1);
        r(end+1) = toc;
    end
    mean(f)/mean(r)
end

end