function markerNames = findMarkerNamesCsv(headerCell)

% All markers in cell delimited by ",,," Assumes all marker names are 4
% chars long so can use matrix instead of cell

% First 2 cols of data are for frame and subframe, then x,y,z coord's for
% each marker in order of names of markers in cell

indSep = strfind(headerCell,',,,');
for i = 1:length(indSep)+1
    if i == length(indSep)+1
        curInd = length(headerCell);
    else
        curInd = indSep(i);
    end
    if i == 1
    	beg = strfind(headerCell(1:curInd),':');
    else
        beg = strfind(headerCell(lastInd:curInd),':') + lastInd - 1;
    end
    lastInd = curInd;
    if ~isempty(beg) && length(headerCell((beg+1):(curInd-1))) == 4
        markerNames(i,:) = headerCell((beg+1):(curInd-1));
    end
end
   