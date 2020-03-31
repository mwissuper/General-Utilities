function dataset = formatSingleTrial(data, delimiters, options, trialStr)
%% FORMATSINGLETRIAL: Reformats the structure of VICON data
%
%   dataset = formatSingleTrial(filename,trialPath) takes a data array
%   DATA containing VICON data, outputted by the processC3DtoMAT function,
%   and converts the structure format.
%
%   FORMATSINGLETRIAL organizes the data by the type of collected data. The
%   fields of the output are:
%           Markers: A structure separating the marker data by subject, and
%           then by marker ID. Any particular marker can be accessed by
%           Marker.(Subject).(MarkerID) so long as you know the name of the
%           marker.
%           PlateMotion: stores the motions of the platform according to
%           their name in platemotionid
%           PlateForce: stores the digital force plate data according to
%           the name in plateforceID
%           EMG: stores EMG according to the EMG box it originated from.
%           Other: contains the data stored in the original 'Other'
%           field
%
%   The goal of the data variable is to make processing easy and quick
%   by consolidating the experimental data into one file with readily
%   accessible components. Much of the work is simply converting IDs into
%   fieldnames so the data can be accessed directly, without the need to
%   referencing the ID array for the location of the data.
%
%   formatSingleTrial can also take two additional arguments, delimiters and
%   options, in the form formatSingleTrial(DATA,delimiters,options).
%       delimiters: when reading in variable names from Vicon data, names
%           can be split into fields and subfields by specifying delimiters
%           as a cell array. The default delimiters are {. :}
%       options: specifying additional options allows concatData to skip
%       over some parts of the data processing. options must be specified
%       as a cell array of strings, with possible arguments:
%           noMarkers, noPlateMotions, noPlateForces, noEMG, noOther,
%           noAnalog, noDigital, noGaitEvents, noAnalogForce
%       arguments for the options are not case sensitive.

%   Luke Drnach
%   March 13, 2019

%% Input Checking
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

if nargin < 4
   trialStr = ''; 
end
% Create an empty structure
dataset = struct();

%% Organize the data

% Process Marker Data
if ~any(strcmpi('nomarkers',options)) && isfield(data,'video') && isempty(data.video.markersid)~=1
    videoFields = fieldnames(data.video);
    dataIdx = strcmpi(videoFields,'markers');
    IDidx = strcmpi(videoFields,'markersID');
    dataset.Markers = Organize(data.video.(videoFields{dataIdx}),data.video.(videoFields{IDidx}),delimiters,options);
    dataset.Markers.samplerate = data.video.samplerate;
else
    disp(['Skipping Marker Data in', trialStr])
end
if isfield(data,'analog') && ~any(strcmpi('noanalog',options))
    analogFields = fieldnames(data.analog);
    %Process PlateMotion Data
    if isfield(data.analog,'platemotion') && isempty(data.analog.platemotion)~=1 && ~any(strcmpi('noplatemotion',options))
        dataIdx = strcmpi(analogFields,'platemotion');
        IDidx = strcmpi(analogFields,'platemotionid');
        dataset.PlateMotion = Organize(data.analog.(analogFields{dataIdx}),data.analog.(analogFields{IDidx}),delimiters,options);
        dataset.PlateMotion.samplerate = data.analog.samplerate;
    else
        disp(['Skipping Plate Motion Data in ' trialStr]);
    end
    % Process EMG Data
    % The EMG Data has little actual structure in it. Will revise this
    % once I have a better way of organizing the EMG data
    if isfield(data.analog,'emg') && isempty(data.analog.emg)~=1 && ~any(strcmpi('noemg',options))
        dataIdx = strcmpi(analogFields,'emg');
        IDidx = strcmpi(analogFields,'emgID');
        dataset.EMG = Organize(data.analog.(analogFields{dataIdx}),data.analog.(analogFields{IDidx}),delimiters,options);
        dataset.EMG.samplerate = data.analog.samplerate;
    else
        disp(['No EMG Data in ',trialStr])
    end
    %Process Analog Force Plate Data
    if isfield(data.analog,'plateforces') && isempty(data.analog.plateforces) ~= 1 && ~any(strcmpi('noAnalogForce',options))
        dataIdx = strcmpi(analogFields,'plateforces');
        IDidx = strcmpi(analogFields,'plateforcesID');
        dataset.AnalogForce = Organize(data.analog.(analogFields{dataIdx}),data.analog.(analogFields{IDidx}),delimiters,options);
        dataset.AnalogForce.samplerate = data.analog.samplerate;
    else
        disp(['Skipping Analog Plate Force Data in ',trialStr]);
    end
    %Process Data from Field 'Other' (Usually force and
    %accelerometry data)
    if isfield(data.analog,'otherid') && isempty(data.analog.otherid)~=1 && ~any(strcmpi('noother',options))
        dataIdx = strcmpi(analogFields,'other');
        IDidx = strcmpi(analogFields,'otherID');
        dataset.Other = Organize(data.analog.(analogFields{dataIdx}),data.analog.(analogFields{IDidx}),delimiters,options);
        dataset.Other.samplerate = data.analog.samplerate;
    else
        disp(['Skipping data labelled Other ',trialStr])
    end
else
    disp(['Skipping analog data in ', trialStr])
end
% Process Gait Event Data
if isfield(data,'GaitEvents') && ~isempty(data.GaitEvents) && ~any(strcmpi('noGaitEvents',options))
    dataset.GaitEvents = data.GaitEvents;
else
    disp(['Skipping Gait Event data in ', trialStr]);
end
%Process Digital Force Plate Force Data
if isfield(data,'digital') && ~any(strcmpi('nodigital',options))
    digitalFields = fieldnames(data.digital);
    if isfield(data.digital,'plateforces') && isempty(data.digital.plateforces)~= 1 && ~any(strcmpi('noplateforces',options))
        dataIdx = strcmpi(digitalFields,'plateforces');
        IDidx = strcmpi(digitalFields,'plateforcesID');
        dataset.PlateForces = Organize(data.digital.(digitalFields{dataIdx}),data.digital.(digitalFields{IDidx}),delimiters,options);
        dataset.PlateForces.samplerate = data.digital.samplerate;
    else
        disp(['Skipping Digital Force Plate Data in ',trialStr])
    end
else
    disp(['Skipping Digital Data in ',trialStr])
end
end

function orgData = Organize(data,dataIDs,delimiters,options)
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
    if isempty(ID) && ~any(strcmpi('noUnknowns', options))
        fname = ['Unknown',num2str(p)];
        orgData.(fname) = sortData;
        p = p+1;
    elseif ~isempty(ID)
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