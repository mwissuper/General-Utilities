function Markers = windowMarkers(Markers,start,stop)
%% WINDOWMARKERS: Returns a subset of VICON Marker data
%
%   MARKERS = windowMarkers(MARKERS,start,stop) returns a structure MARKERS
%   that contains all the fields of the input structure MARKERS. Both the
%   input and output are assumed be structures with array fields. For all
%   nonscalar array fields, windowMarkers returns only the rows of each
%   field between the indices START and STOP. 

%   Luke Drnach
%   December 5, 2018

fields = fieldnames(Markers);
for n = 1:length(fields)
    if ~isscalar(Markers.(fields{n}))
        Markers.(fields{n}) = Markers.(fields{n})(start:stop,:);
    end
end
end