function [dataset] = formatExperimentalData(filename,trialPath,delimiters,options,processedFlag)
%formatExperimentalData: Generates a single working set of data
%
%   [data] = formatExperimentalData(filename,trialPath) reads in an Excel file named
%   filename,containing a list of trial names and other trial information.
%   Trials are read in according to the names in the first column of
%   filename, and are read in from the location trialPath. If trialPath is
%   not specified, trials are read in from the present working directory.
%
%   If filename is a cell array of strings {filename,sheet}, then
%   concatData will read the Excel sheet specified by sheet in the file
%   specified by filename.
%
%   formatExperimentalData can also take two additional arguments, delimiters and
%   options, in the form formatExperimentalData(filename,trialPath,delimiters,options).
%       delimiters: when reading in variable names from Vicon data, names
%           can be split into fields and subfields by specifying delimiters
%           as a cell array. The default delimiters are {. :}
%       options: specifying additional options allows concatData to skip
%       over some parts of the data processing. options must be specified
%       as a cell array of strings, with possible arguments:
%           noMarkers, noPlateMotions, noPlateForces, noEMG, noOther,
%           noAnalog, noDigital, noGaitEvents, noAnalogForce
%       arguments for the options are not case sensitive.
%
%   The data structure array organizes the data by trial - each element
%   of the array is a structure containing that trial's relevant data. The
%   fields of each structure are as follows:
%           Markers: A structure separating the marker data by subject, and
%           then by marker ID. Any particular marker can be accessed by
%           Marker.(Subject).(MarkerID) so long as you know the name of teh
%           marker.
%           PlateMotion: stores the motions of the platform according to
%           their name in platemotionid
%           PlateForce: stores the digital force plate data according to
%           the name in plateforceID
%           EMG: stores EMG according to the EMG box it originated from.
%           Other: contains the data stored in the original 'Other'
%           field
%
%   The structure array also returns a field for each column in the Excel
%   Worksheet, with the data in the field matching the entry in the column
%   in the Worksheet.
%
%   The goal of the data variable is to make processing easy and quick
%   by consolidating the experimental data into one file with readily
%   accessible components. Much of the work is simply converting IDs into
%   fieldnames so the data can be accessed directly, without the need to
%   referencing the ID array for the location of the data.

%   Luke Drnach
%   November 12, 2016
%   Updated: November 19, 2016
%           - Updated handling of force handle to generalize to 'other'
%             field
%           - Included fields for EMG data
%           - Included fields for acceleromery data
%   Updated January 25, 2017
%           - Excel file and trial data need not be on the same path
%           - Extra handling for data - elimination of all characters that
%           cannot be used as field names
%           - Import options now available: noMarkers, noEMG, noPlateMotion,
%           noPlateForces, noOther
%           - Added TrialInfo field. Any data stored in the read-in
%           worksheet is imported into TrialInfo
%   Updated February 24, 2017
%           - Can import sheets of an excel file by specifying filename as
%           a cell array, where the first element is the file and the
%           second is the sheet.
%   Updated April 3, 2017
%           - Now checks for GaitEvents, and copies over GaitEvents
%           structures into the final output. Can specify 'noGaitEvents'
%           under options to skip GaitEvent processing.
%   Updated April 12, 2017
%           -Now checks for analog force plate readings, and copies analog
%           forces into AnalogForce in the output. Specifying
%           'noAnalogForce' under options skips this step.
%   Updated June 7, 2017
%           -Final fieldnames in original data structure are now
%           case-insensitive. For example, concatData can now process both
%           data of the format data.video.markersID or data.video.markersid.
%           Similar functionality holds for all fields.
%   Updated October 26, 2017
%           -Specifying option 'noUnknowns' prevents concatData from adding
%           fields with no valid name to the output as 'Unknownx', where x
%           is a number.
%   Updated March 13, 2018
%           -Refactored to make a separate function for organizing
%           individual trials. 
%           -Renamed to 'formatExperimentalData'
%           -Allows users to select either the 'rawData' or the processed
%           'data' from VICON. Default is 'rawData'

% INPUT CHECKING
narginchk(1,5);
flag = 0;
if isempty(filename) || ~ischar(filename)
    if iscell(filename)
        flag = 1;
    else
        error('Filename must be either a string or a cell array of strings')
    end
end
if nargin > 1
    if isempty(trialPath)
        trialPath = pwd;
    elseif ~ischar(trialPath)
        error('trialPath must be either a character array or empty')
    end
    if nargin > 2
        if isempty(delimiters)
            delimiters = {':','.'};
        elseif ~iscell(delimiters)
            error('Delimiters must be a cell array');
        end
        if nargin > 3
            if isempty(options)
                options = {''};
            elseif ~iscell(options)
                error('options must be a cell array of strings');
            end
        else
            options = {''};
        end
    else
        delimiters = {'.',':'};
        options = {''};
    end
else
    trialPath = pwd;
    delimiters = {'.',':'};
    options = {''};
end
if nargin < 5 || isempty(processedFlag)
    processedFlag = 0;
else
    processedFlag = 1;
end
%Check that the excel file exists in the present directory
if flag
    if exist(filename{1},'file') == 2
        [~,~,raw] = xlsread(filename{1},filename{2});
        DataTable = cell2table(raw(2:end,:),'VariableNames',raw(1,:));
    else
        error('File does not exist');
    end
else
    
    if exist(filename,'file') == 2
        DataTable = readtable(filename);
    else
        error('File does not exist')
    end
end
%Initialize the data structure
[NTrials,NInfo] = size(DataTable);
dataset(1:NTrials) = struct('Info',[]);
emptysets = false(1,NTrials);  %Keep a record of which trials are not included in the working set
VarNames = DataTable.Properties.VariableNames;
%Read in individual trials and store them in the data variable
for n = 1:NTrials
    %Load the files listed in the first column of DataTable
    trialStr = DataTable{n,1}{:};
    if exist([trialPath,'\',trialStr,'.mat'],'file') ~= 2
        %Checks if the file exists; if it does not the row of the structure
        %is empty and the row number is recorded so it can be removed later
        disp([trialStr,' is not in ',trialPath]);
        emptysets(n) = true;
    else
        %If the trial exists, the trial is loaded and the information in
        %DataTable is recorded in Info according to the column headings.
        %Then the rest of the trial data is recorded.
        disp(['Processing ',trialStr]);
        if processedFlag
            load([trialPath,'\',trialStr,'.mat'],'data');
        else
            load([trialPath,'\',trialStr,'.mat'],'rawData');
            data = rawData;
        end
        for p = 1:NInfo
            if iscell(DataTable{n,p})
                dataset(n).Info.(VarNames{p})=DataTable{n,p}{:};
            else
                dataset(n).Info.(VarNames{p})=DataTable{n,p};
            end
        end
        % Format the individual trial
        formatData = formatSingleTrial(data, delimiters, options, trialStr);
        % Copy the formatted trial fields into the main experimental
        % dataset
        fields = fieldnames(formatData);
        for k = 1:length(fields)
            dataset(n).(fields{k}) = formatData.(fields{k});
        end
    end
end
% Remove slots for trials that do not exist in the workspace
dataset(emptysets) = [];
end