function [dataset] = concatDataMat(filename,trialPath,delimiters,options)
%concatData: Generates a single working set of data 
%
%   [data] = concatData(filename,trialPath) reads in a mat file named 
%   filename,containing a list of trial names and other trial information. 
%   Trials are read in according to the names in the first column of
%   filename, and are read in from the location trialPath. If trialPath is
%   not specified, trials are read in from the present working directory.
%
%   If filename is a cell array of strings {filename,sheet}, then
%   concatData will read the Excel sheet specified by sheet in the file
%   specified by filename.
%
%   concatData can also take two additional arguments, delimiters and
%   options, in the form concatData(filename,trialPath,delimiters,options).
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
%           TrialInfo: Trial information contained in the Excel Worksheet
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
%           Every field except TrialInfo also contains a variable
%           samplerate giving the sampling rate (Hz) of the data in that
%           field.
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

% INPUT CHECKING
narginchk(1,4);
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
        temp = load(filename);
        DataTable = temp.data;
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
       load([trialPath,'\',trialStr,'.mat']);
       for p = 1:NInfo
           if iscell(DataTable{n,p})
               dataset(n).Info.(VarNames{p})=DataTable{n,p}{:};
           else
               dataset(n).Info.(VarNames{p})=DataTable{n,p};
           end
       end
       % Process Marker Data
       if isfield(data,'video') && isempty(data.video.markersid)~=1 && ~any(strcmpi('nomarkers',options))         
           dataset(n).Markers = Organize(data.video.markers,data.video.markersid,delimiters);
           dataset(n).Markers.samplerate = data.video.samplerate;
       else
           disp(['Skipping Marker Data', trialStr])
       end
       if isfield(data,'analog') && ~any(strcmpi('noanalog',options))
           %Process PlateMotion Data
           if isfield(data.analog,'platemotion') && isempty(data.analog.platemotion)~=1 && ~any(strcmpi('noplatemotion',options))
               dataset(n).PlateMotion = Organize(data.analog.platemotion,data.analog.platemotionid,delimiters);
               dataset(n).PlateMotion.samplerate = data.analog.samplerate;
           else
               disp(['Skipping Plate Motion Data in ' trialStr]);
           end
           % Process EMG Data
           % The EMG Data has little actual structure in it. Will revise this
           % once I have a better way of organizing the EMG data
           if isfield(data.analog,'emg') && isempty(data.analog.emg)~=1 && ~any(strcmpi('noemg',options))
               dataset(n).EMG = Organize(data.analog.emg,data.analog.emgid,delimiters);
               dataset(n).EMG.samplerate = data.analog.samplerate;
           else
               disp(['No EMG Data in ',trialStr])
           end
           %Process Analog Force Plate Data
           if isfield(data.analog,'plateforces') && isempty(data.analog.plateforces) ~= 1 && ~any(strcmpi('noAnalogForce',options))
               dataset(n).AnalogForce  = Organize(data.analog.plateforces,data.analog.plateforcesid,delimiters);
               dataset(n).AnalogForce.samplerate = data.analog.samplerate;
           else
               disp(['Skipping Analog Plate Force Data in ',trialStr]);
           end
           %Process Data from Field 'Other' (Usually force and
           %accelerometry data)
           if isfield(data.analog,'otherid') && isempty(data.analog.otherid)~=1 && ~any(strcmpi('noother',options))
               dataset(n).Other = Organize(data.analog.other,data.analog.otherid,delimiters);
               dataset(n).Other.samplerate = data.analog.samplerate;
           else
               disp(['Skipping data labelled Other ',trialStr])
           end
       else
           disp(['Skipping analog data in ', trialStr])
       end
       % Process Gait Event Data
       if isfield(data,'GaitEvents') && ~isempty(data.GaitEvents) && ~any(strcmpi('noGaitEvents',options))
           dataset(n).GaitEvents = data.GaitEvents;
       else
           disp(['Skipping Gait Event data in ', trialStr]);
       end
       %Process Digital Force Plate Force Data
       if isfield(data,'digital') && ~any(strcmpi('nodigital',options))
           if isfield(data.digital,'plateforces') && isempty(data.digital.plateforces)~= 1 && ~any(strcmpi('noplateforces',options))
              dataset(n).PlateForces = Organize(data.digital.plateforces,data.digital.plateforcesid,delimiters);
              dataset(n).PlateForces.samplerate = data.digital.samplerate;
           else
               disp(['Skipping Digital Force Plate Data in ',trialStr])
           end
       else
           disp(['Skipping Digital Data in ',trialStr])
       end     
  end
end
dataset(emptysets) = [];
end

function orgData = Organize(data,dataIDs,delimiters)
charAllow = 'qwertyuiopasdfghjklzxcvbnm1234567890_';    %Characters allowed in field names
charAllow = [charAllow,strjoin(delimiters,'')];         %Includes delimiters to make searching easier.
orgData = struct();
% Convert to cell array: Removes leading and trailing whitespace 
if ~iscell(dataIDs)
    dataIDs = cellstr(dataIDs);
end
K = length(dataIDs);
p = 1;
N = ndims(data);
%Define two anonymous functions to make sure delimited field names are
%valid
fun1 = @(X) isletter(X(1));
fun2 = @(X) horzcat('a',X);
for k=1:K
    %Pull out the data based on the input dimensions
    if N == 2
        sortData = data(:,k);
    else
        sortData = squeeze(data(:,k,:));
    end
    %Get the ID tag as a character array
    ID = dataIDs{k};
    %Check that the ID starts with a letter
    if ~isletter(ID(1))
        idx = isletter(ID);
        idx1 = find(idx==1,1,'first');
        ID = ID(idx1:end);
    end
    %Now remove everything but the allowable characters
    ID(~ismember(lower(ID),charAllow)) = [];
    if isempty(ID)
        fname = ['Unknown',num2str(p)];
        orgData.(fname) = sortData;
        p = p+1;
    else
       %Divide the ID tag based on the delimiters
       fnames = strsplit(ID,delimiters);
       %Check that all delimited field names are valid field names
       chk = cellfun(fun1,fnames);
       if any(chk)
          fnames(~chk) = cellfun(fun2,fnames(~chk),'UniformOutput',0);  
       end
       %Store the data in the final structure
       orgData = setfield(orgData,fnames{:},sortData);
    end  
end
end 