function [TrialData] = processHHI(TrialData)
%PROCESSHHI: Data processing for HHI experiment
%
%   TrialData = processHHI(TrialData) performs a series of analyses on
%   on the data stored in workset. processHHI loops through the data and
%   processes trials of type 'Assist Ground' or 'Assist Beam'. The analyses
%   currently performed include:
%       realigning the forces to match the Vicon coordinate system
%       WorkLoop analysis to estimate work removed during the trial
%   
%   processHHI adds a field Results to TrialData which contains:
%       Forces: the realigned interaction forces from the trial (N)
%       AssistArm: the arm vector for the assistance provider (mm)
%       BeamArm:   the arm vector for the beam walker (mm)
%       AssistWork: the work removed from the assistance provider (mJ)
%       BeamWork:   the work removed from the beam walker (mJ)
%       time:    the time vector for the forces, arm vectors, and work
%       startidx: the collection index for which analysis starts
%       stopidx: the collection index for which analysis ends
%
%   NOTE: much of the current data is yet to be converted to the proper
%   unit system and calibrated. Attempts to read results in any given unit
%   system would be ill-advised.

% Luke Drnach
% May 30, 2017
% Neuromechanics Lab
% Georgia Tech and Emory University
    
    %% Initialization
    % Number of trials in the dataset
    NTrials = length(TrialData); %Number of trials in the dataset
    % Load the force handle standard for realigning the forces
    Ref = genHandleRef('Stat03.mat');
    % Force gain: Convert voltage readings to Newtons
    Fgain = 105.6;
    %% Trial Analysis Main Loop
    for n = 1:NTrials
       disp(['Processing ',TrialData(n).Info.Trial]);
        % Work Analysis (Assisted Trials ONLY)
       if any(strcmp(TrialData(n).Info.Condition,{'Assist Beam','Assist Ground'}))
           % Now compute the start and stop time of the trial.
           % Start time is computed from the vertical displacement of the Person
           %on the Beam's Left Heel.
           %Start time is determined by the location of the first half-maximum
           %of the vertical displacement of the heel. The half-maximum is
           %computed by first finding all maxima with prominence above half the
           %range of the vertical displacement, and then finding the location of
           %the first point where the displacement is greater than the
           %half-maximum.
           LHeel = medfilt1(TrialData(n).Markers.PersononBeam.LHEE(:,3));
           peak = findpeaks(LHeel,'MinPeakProminence',range(LHeel)/2);
           startidx = find(LHeel >= peak(1)/2,1,'first');
           time = (0:length(LHeel)-1)./TrialData(n).Markers.samplerate;
           % Get the Shoulders and Fingers for both subjects
           AssistRSHO = TrialData(n).Markers.AssistanceProvider.RSHO;
           AssistRFIN = TrialData(n).Markers.AssistanceProvider.RFIN;
           BeamerLSHO = TrialData(n).Markers.PersononBeam.LSHO;
           BeamerLFIN = TrialData(n).Markers.PersononBeam.LFIN;
           % Check that all markers are present and adjust start time
           % accordingly
           if any(abs(AssistRSHO(:,3)) < 50)
               idx = find(abs(AssistRSHO)>50,1,'first') + 20;
               startidx = max([startidx,idx]);
           end
           if any(abs(AssistRFIN(:,3)) < 50)
               idx = find(abs(AssistRFIN)>50,1,'first') + 20;
               startidx = max([startidx,idx]);
           end
           if any(abs(BeamerLSHO(:,3)) < 50)
               idx = find(abs(BeamerLSHO)>50,1,'first') + 20;
               startidx = max([startidx,idx]);
           end
           if any(abs(BeamerLFIN(:,3)) < 50)
               idx = find(abs(BeamerLSHO)>50,1,'first') + 20;
               startidx = max([startidx,idx]);
           end
           %End time is computed as the time when the CLAV marker has travelled
           %12 ft (length of the beam) for all assisted trials.
           ClavY = TrialData(n).Markers.PersononBeam.CLAV(startidx:end,2);
           ClavY = (ClavY-ClavY(1))./304.8; %Convert mm to ft
           stopidx = find(ClavY >=12,1,'first');
           if isempty(stopidx)
               stopidx = length(ClavY) + startidx - 1;
           end
           % Compute the arm vectors from the trial data. The arm is  
           % defined as the vector pointing from the shoulder of the
           % participant to the hand. 
           % Median filtering is used to remove large single point jumps in
           % the data
           AssistArm = medfilt1(AssistRFIN,[],[],1) - medfilt1(AssistRSHO,[],[],1);
           BeamArm = medfilt1(BeamerLFIN,[],[],1) - medfilt1(BeamerLSHO,[],[],1);          
           %Now use only the parts in the defined trial time windows
           AssistArm = AssistArm(startidx:stopidx,:);
           BeamArm = BeamArm(startidx:stopidx,:);                    
           
           %Now get the handle forces and resample them to match the marker
           %sampling rate 
           Force = zeros(3,length(time));
           Field = {'Fx','Fy','Fz'};
           MarkerRate = TrialData(n).Markers.samplerate;
           ForceRate = TrialData(n).Other.samplerate;
           for k = 1:3
                Force(k,:) = resample(TrialData(n).Other.ElectricPotential.(Field{k}),MarkerRate,ForceRate);  
           end
           %Median filter forces to reduce noise
           Force = medfilt1(Force,[],[],2);
           %Realign forces with the Vicon Frame
           Force = alignForce(Force,TrialData(n).Markers.ForceHandle,Ref);
           %Apply gain and limit forces to trial time window
           Force = Fgain*Force(:,startidx:stopidx);
           time = time(startidx:stopidx)-time(startidx);
           
           %%%%% WORK ANALYSIS %%%%%
           % First, transpose the Arm vectors so they are columns
           AssistArm = AssistArm';
           BeamArm = BeamArm';
           % Next, compute the change in arm length between time points
           AssistDiff = diff(AssistArm,1,2)*MarkerRate;
           BeamDiff = diff(BeamArm,1,2)*MarkerRate;
           % Now, calculate Power as the dot product of force
           % and displacement (the 1000 is to convert mW to W)
           AssistPow = [0,dot(-Force(:,2:end),AssistDiff)]./1000;
           BeamPow = [0,dot(Force(:,2:end),BeamDiff)]./1000;
           AssistWork = cumsum(AssistPow)./MarkerRate;
           BeamWork = cumsum(BeamPow)./MarkerRate;          
           % Store Results
           TrialData(n).Results.AssistArm = AssistArm;
           TrialData(n).Results.BeamArm = BeamArm;
           TrialData(n).Results.Force = Force;
           TrialData(n).Results.AssistWork = AssistWork;
           TrialData(n).Results.BeamWork = BeamWork;
           TrialData(n).Results.AssistPow = AssistPow;
           TrialData(n).Results.BeamPow = BeamPow;
           TrialData(n).Results.Clav = TrialData(n).Markers.PersononBeam.CLAV(startidx:stopidx);
           TrialData(n).Results.time = time;
           TrialData(n).Results.startidx = startidx;
           TrialData(n).Results.stopidx = stopidx;
       end
    end
end

function [Ref] = genHandleRef(Static)
%Generates a Reference for aligning the forces from a static file
S = load(Static);       %Load the static file
Markers = S.Markers;
IDs = S.MarkerID;
[~,Idx] = max(size(Markers));
Markers = squeeze(mean(Markers,Idx));
[row,col] = size(Markers);
if row > col
    Markers = Markers';
end
Ref = zeros(3,4);
for n = 1:4
    ID = IDs(n,:);
    ID(ismember(ID,' ')) = [];
    %Organize the static markers where the row is (xyz) and the cols are
    %  (frontleft, frontright, backleft, backright)
    if strcmp(ID,'frontleft')
        Ref(:,1) = Markers(:,n);
    elseif strcmp(ID,'frontright')
        Ref(:,2) = Markers(:,n);
    elseif strcmp(ID,'backleft')
        Ref(:,3) = Markers(:,n);
    else
        Ref(:,4) = Markers(:,n);
    end
end
%Translate all rows such that the front right is origin
Ref = Ref - repmat(Ref(:,1),1,4);
Ref = Ref(:,2:end);     %Remove the frontright / origin
%Generate a rotation which aligns the front left with the x-axis WITHOUT
%rotation about the x-axis
theta = atan(Ref(3,1)./Ref(1,1));
R1  = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Ref = R1*Ref;   %Rotates about absolute y-axis to eliminate absolute z-component
phi = atan(Ref(2,1)./Ref(1,1));
R2 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
%Final rotation about absolute z-axis brings front left into alignment
Ref = R2*Ref;
%If necessary, reflect about the yz plane
if Ref(1,1) < 0
   R3 = eye(3);
   R3(1,1) = -1;
   Ref = R3*Ref;
end
%The assumption is that, in this orientation, the force sensor is aligned
%with the global (Vicon) coordinates. The problem remaining is to bring all
%other force marker sets into alignment with the reference.
end

function [Force] = alignForce(Force,Markers,Ref)

N = length(Force);
%For each frame, generate the rotation that brings the handle into
%alignment with the reference. Then apply that rotation to the force for
%that frame.
for n = 1:N
    %Get the Handle Markers for this frame
    Handle = [Markers.frontleft(n,:)',Markers.frontright(n,:)',Markers.backleft(n,:)',Markers.backright(n,:)'];
    Rot = alignHandle(Handle,Ref);  %Generate the rotation matrix
    Force(:,n) = Rot*Force(:,n);    %Rotate the forces.
end
end

function [Rot] = alignHandle(Handle,Ref)
%Handle and Ref must both be stored such with the rows corresponding to
%(x,y,z) and the cols corresponding to the markers (front left, front
%right, back left, back right)

%Step 1: Translate the force handle s.t. front right is at origin
Handle = Handle - repmat(Handle(:,1),1,4);
Handle = Handle(:,2:end);     %Remove the origin from the matrix
%Step 2: Align the front left with the x-axis
R = genVecAlign(Handle(:,1),[1 0 0]');
Handle = R*Handle;
%Step 3: Align the Normals of the plane defined by the force handle markers
HandleNorm = cross(Handle(:,1),Handle(:,2));
RefNorm = cross(Ref(:,1),Ref(:,2));
R2 = genVecAlign(HandleNorm,RefNorm);
Rot = R2*R;
end

function R = genVecAlign(x,y)
%Normalize X,Y, so the resulting rotation matrix is orthonormal
x = x./norm(x);     %Normalize x,y
y = y./norm(y);
v = cross(x,y);     %Cross product
c = dot(x,y);       %Cosine of the angle between x,y
if c ~= -1          %If the vectors are not perfectly opposite one another
    V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + V + V^2*(1/(1+c));
else                %If the vectors are complete opposites, the rotation is a coordinate inversion (reflection)
    R = -eye(3);
end
end
