% Find lag as delay that results in max xcorr. Do not allow anticorr as
% data is cyclic

function [m, curlag] = xcorrLag(x,y,maxlag)

 % Must remove DC before do xcorr
 xdc = x - nanmean(x);
 ydc = y - nanmean(y); 
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