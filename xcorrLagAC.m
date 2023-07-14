% Find lag as delay that results in max xcorr. Allow anticorr 

function [m, curlag] = xcorrLagAC(x,y,maxlag)

 % Must remove DC before do xcorr
 xdc = x - nanmean(x);
 ydc = y - nanmean(y); 
 if isnan(maxlag)
     [r,lags] = xcorr(xdc,ydc,'normalized'); 
 else
     [r,lags] = xcorr(xdc,ydc,maxlag,'normalized'); 
 end

 [temp,imax] = max(abs(r));
 curlag = lags(imax);
 
 if r(imax) < 0
     m = - temp;
 else
     m = temp;
 end