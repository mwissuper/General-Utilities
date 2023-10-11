% Find lag as delay that results in max xcorr. Do not allow anticorr as
% data is cyclic. curlag is in samples

function [m, curlag] = xcorrLag(x,y,maxlag)

 % Must remove nan's before do xcorr
 ix = find(~isnan(x)); iy = find(~isnan(y));
 % Must remove DC before do xcorr
 ind = intersect(ix,iy);
 xdc = x(ind) - nanmean(x);
 ydc = y(ind) - nanmean(y); 
 if isnan(maxlag)
     [r,lags] = xcorr(xdc,ydc,'normalized'); 
 else
     [r,lags] = xcorr(xdc,ydc,maxlag,'normalized'); 
 end

 [m,imax] = max(r);
 curlag = lags(imax);

 if m < 0 % anti-correlation, can't use data
    m = nan;
    curlag = nan;
 end