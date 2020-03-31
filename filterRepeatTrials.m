function trials = filterRepeatTrials(trials)
%% FILTERREPEATTRIALS: Removes repeated trials from the data
%
%   filtTrials = filterRepeatTrials(trials) searches the input TRIALS for
%   repeated trials, indicated by a NOTES field, and replaces the original
%   trial with its repeated trial. The repeated trial is then removed from
%   the TRIALS data structure, and returned as FILTTRIALS.
%
%   TRIALS should be a table of trial statistics (from getPointStatsHHI).

%   Luke Drnach
%   December 7, 2018


if isa(trials,'table')
    % Pull the NOTES for all the trials
    notes = trials.Notes;
    % Search for NOTES that start with 'REPEAT'
    filtIdx = find(strncmp(notes,'Repeat',6));
    if isempty(filtIdx)
        % If there are no repeats found, return the input data structure.
        return;
    else
        replaced = true(1,length(filtIdx));
        remIdx = zeros(1,length(filtIdx));
        % If there are repeats, we loop over the number of repeats
        for n = 1:length(filtIdx)
            % Parse the trial number from the repeat string
            repStr = trials.Notes{filtIdx(n)};       %String containing 'Repeat'
            numIdx = strfind(repStr,'Trial') + 5;    %Somewhere there should be 'Trial##'. Find the location of the first #.
            tNum = str2double(repStr(numIdx:end));   %Trial Number to replace
            % Get the associated subject number
            subNum = trials.Subject(filtIdx(n));
            % Find the matching subject and trial
            repIdx = find(and(trials.Subject == subNum, trials.TrialNumber == tNum));
            % Copy the repeated trial into the place of the original trial
            if ~isempty(repIdx)
                remIdx(n) = repIdx;
            else
                replaced(n) = false;
            end
        end
        % Delete all the repeat trials
        trials(remIdx(replaced), :) = [];
    end
else
    Info = [trials.Info];
    notes = {Info.Notes};
    trialNum = {Info.Trial};
    subject = {Info.Subject_1};
    % Convert the trial numbers from a cell array into a numeric array
    trialNumber = nan(1,length(trialNum));
    for n = 1:length(trialNum)
       if strncmp(trialNum{n},'Trial',5)
          trialNumber(n) = str2double(trialNum{n}(6:end)); 
       end
    end
    % Search for NOTES that start with 'REPEAT'
    filtIdx = find(strncmp(notes,'Repeat',6));
    % Replace Trials that start with STAT with NAN
    if isempty(filtIdx)
       return; 
    else
        replaced = true(1,length(filtIdx));
        remIdx = zeros(1,length(filtIdx));
        for n = 1:length(filtIdx)
            % Get the Subject for this trial
            repSub = trials(filtIdx(n)).Info.Subject_1;
            % Parse the trial number from the repeat string
            repStr = notes{filtIdx(n)};              %String containing 'Repeat'
            numIdx = strfind(repStr,'Trial') + 5;    %Somewhere there should be 'Trial##'. Find the location of the first #.
            tNum = str2double(repStr(numIdx:end));   %Trial Number to replace
            % Find the matching subject and trial
            repIdx = find(and(strcmp(repSub,subject), trialNumber == tNum));
            % Copy the repeated trial into the place of the original trial
            if ~isempty(repIdx)
                remIdx(n) = repIdx;
            else
                replaced(n) = false;
            end
        end
    end
    % Delete all repeat trials
    trials(remIdx(replaced)) = [];
end
end