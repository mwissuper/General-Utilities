function [TrialData] = relabelData(TrialData,field,oldnames,newnames)
    %RELABELDATA: relabels fields within the TrialData data structure
    %
    %   TrialData = relabelData(TrialData,field,oldnames,newnames) converts
    %   TrialData(n).field.oldnames into
    %   TrialData(n).field.newnames. oldnames and newnames can be cell
    %   arrays specifying multiple fields; however, the lengths of oldnames
    %   and newnames must match. FIELD is a single string and cannot be
    %   a cell array. For nonscalar structure TrialData, RELABELDATA
    %   operates on every element in the structure array.   
    
    %   Luke Drnach
    %   May 31, 2017
    %   Neuromechanics Lab
    %   Georgia Tech and Emory
    
    %% Input Check and Initialization
    if ~iscell(oldnames)
       oldnames = {oldnames}; 
    end
    if ~iscell(newnames)
       newnames = {newnames}; 
    end
    if ~iscell(field)
        field = {field};
    end
    NOld = length(oldnames);
    NNew = length(newnames);
    if NOld ~= NNew
       error('There must be as many old field names as new ones'); 
    end
    NTrials = length(TrialData);
    % Set-up the necessary structure S for subsasgn
    S = struct('type',[],'subs',[]);
    for n = 1:length(field)
        S(n).type = '.';
        S(n).subs = field{n};
    end
    %% Relabelling loop
    if ~isfield(TrialData,field)
       disp([field,' is not a field of TrialData']);
       return
    end
    for n = 1:NTrials
          try
              Temp = getfield(TrialData,{n},field{:},{':'});
              for m = 1:NOld
                  if isfield(Temp,oldnames{m})
                    Temp.(newnames{m}) = Temp.(oldnames{m});
                    Temp = rmfield(Temp,oldnames{m});
                  end
              end
              TrialData(n) = subsasgn(TrialData(n),S,Temp);
          catch ME
              %throw(ME);
          end
           %if isfield(TrialData(n).(field),oldnames{m})
              %To relabel, copy the old data into the new field, then
              %delete the old field.
             %TrialData(n).(field).(newnames{m}) = TrialData(n).(field).(oldnames{m});
             %TrialData(n).(field) = rmfield(TrialData(n).(field),oldnames{m});
          %end
    end
end

