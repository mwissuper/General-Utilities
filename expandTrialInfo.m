function [newTrialData] = expandTrialInfo(TrialData)
%% EXPANDTRIALINFO: Copies data in TrialInfo field into separate subfields
%
%   TrialData = expandTrialInfo(TrialData) copies the information in the
%   Info field of TrialData, and experimental record created by
%   formatExperimentalData, into separate fields of TrialData. Info is
%   then removed from the TrialData structure. The fields of TrialData that
%   are not Info are unaltered.

if isfield(TrialData,'Info')
    newTrialData = struct();
    for n = 1:length(TrialData)
        % Expand the TrialInfo Field
        infoFields = fieldnames(TrialData(n).Info);
        for k = 1:length(infoFields)
           newTrialData(n).(infoFields{k}) =  TrialData(n).Info.(infoFields{k});
        end
        % Copy the remaining fields
        fields = fieldnames(TrialData(n));
        infoIdx = strcmpi(fields,'Info');
        fields(infoIdx) = [];
        for k = 1:length(fields)
            newTrialData(n).(fields{k}) = TrialData(n).(fields{k});
        end 
    end
else
    % If there's no INFO field, return the input
    newTrialData = TrialData;
end
end

