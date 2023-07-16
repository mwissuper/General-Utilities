function x_filt = filtIgnoreNan(B,A,order,x)

% find each nan gap, and filter data as sections around that gap. Then
% re-insert the nan's after filtering. Butterworth lowpass. Filtfilt
% accounts for transients at beginning and end of data segment, so not
% going to worry about reflection and extrapolation stuff.

% Find all beginnings of gaps
offset = x(2:end);
indNan = isnan(x);
indStop = [];
if isnan(x(1)) == 1 % If start with gap
    ind2(1) = 1;
    temp = find(indNan(1:end-1)==0 & isnan(offset)==1); % index before first nan element of gap
    ind2(2:length(temp)+1) = temp;
else
    ind2 = find(indNan(1:end-1)==0 & isnan(offset)==1); % index before first nan element of gap
end

x_filt = x;
ind1 = 1;
for i = 1:length(ind2) % for number of nan gaps
%     ind1 = ind1
%     temp = ind2(i)
    if length(ind1:ind2(i)) > 3*order % filtfilt only works if this is true
        x_filt(ind1:ind2(i)) = filtfilt(B,A,x(ind1:ind2(i)));
    else
        msg = sprintf('gap too short for filter order %i',order);
%         disp(msg);
    end
    indStop = find(indNan(ind2(i)+1:end)==0,1,'first')+ind2(i)-1; % find end of this gap, the last NaN value in gap. 
    ind1 = indStop+1; % to use for next loop
end
if isempty(indStop) == 0 && length(x_filt(ind1:end)) > 3*order % if there is an end of the last gap (and enough pts in data), need to filter rest of non-nan values after last gap
    x_filt(ind1:end) = filtfilt(B,A,x(ind1:end));
end
