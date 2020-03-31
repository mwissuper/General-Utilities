function [TrialData] = getTrials(TrialData,Condition,Type)
%GETTRIALS: Returns a structure of trials of a particular subtype
%
%   TrialData = getTrials(TrialData,Condition,Type) returns a structure of
%   trials whose experimental condition matches those specified by the pair
%   Condition, Type. For example:
%
%       TrialData = getTrials(TrialData,'Condition','PWS') returns a
%       structure of trials whose Condition is labeled as PWS.
%   
%   In general, getTrials returns all data for which Type =
%   TrialData.Info.(Condition)
%
%   Requirements: TrialData must have a field called Info, which is a
%   structure containing information about an individual trial. Condition
%   must be a subfield of Info, otherwise getTrials returns an empty
%   structure.
%
%
%   Luke Drnach
%   March 14, 2017
%   Georgia Institute of Technology

%   Check that Info is a subfield of TrialData
if ~isfield(TrialData,'Info')
   error('Data structure must have an Info field'); 
end
%   Pull out the trial experimental conditions
Info = [TrialData(:).Info];
if ~isfield(Info,Condition)
   error(['No condition ', Condition]); 
end
ExpConds = {Info(:).(Condition)};
%Compare the experimental conditions against the selected type
idx = ismember(ExpConds,Type);
TrialData = TrialData(idx);
end

