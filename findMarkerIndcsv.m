function markerName = findMarkerIndcsv(headerCell)

% Take input string s of name of marker and find which index it is in
% header.

% Separate all markers from string in cell using ",,,"

% First 2 cols of data are for frame and subframe, then x,y,z coord's for
% each marker in order of names of markers in cell
indSep = strfind(headerCell,',,,');
for i = 1:length(indSep)+1
    if i == length(indSep)+1
        curInd = length(indSep);
    else
        curInd = indSep(i);
    end
    if i == 1
    	beg = strfind(headerCell(1:curInd),':');
        markerName{i} = headerCell((beg+1):(curInd-1));
    else
        beg = strfind(headerCell(lastInd:curInd),':') + lastInd - 1;
    end
    lastInd = curInd;
    markerName{i} = headerCell((beg+1):(curInd-1));
end
   