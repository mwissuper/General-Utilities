% Find lag as delay that results in max xcorr. Allow anticorr 

function [m, curlag] = xcorrLagAC(x,y,maxlag)

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

 [temp,imax] = max(abs(r));
 curlag = lags(imax);
 
 if r(imax) < 0
     m = - temp;
 else
     m = temp;
 end