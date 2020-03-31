function Markers = transposeMarkers(Markers)
%% TRANSPOSEMARKERS: Transposes all the fields of a structure
%
%   Markers = transposeMarkers(Markers) takes a structure MARKERS whose
%   fields are arrays and returns the structure with all the fields
%   transposed. For a MARKERS structure containing VICON marker data,
%   TRANSPOSEMARKERS converts all fields between row and column vector
%   formats.

%   Luke Drnach
%   December 4, 2018

fields = fieldnames(Markers);
for n = 1:length(fields)
    Markers.(fields{n}) = Markers.(fields{n})';
end
end