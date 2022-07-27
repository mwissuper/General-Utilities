function [TrialData] = mainWorkPowerAnalysisMW(inputData,subj,plotCheck)
    %MAINWORKPOWERANALYSIS: Processing for HHI Experimental Data collected after
    %August 2017
    %
    %   TrialData = mainWorkPowerAnalysis(inputData) performs a series of analyses on
    %   the data stored in inputData. mainWorkPowerAnalysis loops through the data
    %   stored in inputData and processes trials of type 'Assist
    %   Ground', 'Assist Beam', or 'Assist Solo.' The analysis includes:
    %       1) Recovering the interaction forces in Vicon Coordinates
    %       2) Calculating power and related metrics
    %   mainWorkPowerAnalysis also processes trials of type 'Solo Beam'; however,
    %   these trials do not require a work analysis
    %
    %   mainWorkPowerAnalysis copies inputData to TrialData and adds a field called Results to TrialData. 
    %   Note: Force calibration and conversion is performed in a separate
    %   file recoverForces.m. Adjustments to recoverForces.m may be needed
    %   as the force sensor is recalibrated. Additionally, recoverForces.m
    %   does not take into account the weight of the force handle, which
    %   may be of imporance later on.
    
    %   Luke Drnach
    %   October 30, 2017
    
    %   Refactored on December 4, 2018. Individual routines were separated
    %   into functions for code clarity and potential future re-use.
    
    % Significantly altered by MW 2019-2020 to perform more analyses and
    % change the way the analysis window is found for beam-walking trials
    
    %% Initialization
    TrialData = inputData; % Copy over all orig data and then can remove bad trials from this output struct. Do all further analysis using data in output struct to keep indices correct.

    % Trial Types to Process
    to_process = {'Assist Beam','Assist beam','Assist Ground','Solo Beam','Solo Ground','Assist Solo'};
    % Assist Trial Names
    assisted = {'Assist Beam','Assist beam','Assist Ground','Assist Solo'};
    % Field names for force and torque
    force_names = {'Fx','Fy','Fz'};
    torque_names = {'Mx','My','Mz'};
    % Butterworth 3rd order lowpass filter, cutoff 10Hz
    [z, p, k] = butter(3,10/50);
    [sos, g] = zp2sos(z, p, k);
    
    % Mass of FT sensor
    mFT = 4/9.81; % (kg)
    
    % For check plots of Assist Beam and Solo Beam only
    if subj == 12 && (plotCheck == 4 || plotCheck == 6)
        numrows = 5;
    else
        numrows = 4; 
    end
    numcols = 5;
    plotind = 0; plot2ind = 0; % 2 figures, one for overground walking and one for beam-walking
    
    if subj == 3
        TrialData(30) = []; % bad processing in Nexus
    elseif subj == 4 
        TrialData(16) = []; % bad processing in Nexus
    elseif subj == 5 
        TrialData(4) = []; % didn't record beginning of trial, must remove later trial before earlier or indices messed up
        TrialData(2) = []; % huge outlier, leaned a lot on partner
    elseif subj == 7 
        TrialData(51) = []; % trial to repeat trial 23 bc no gait mat in trial 23. But doesn't matter for analysis, trial 23 looks fine.
    elseif subj == 14
        TrialData(26) = []; % notes say to skip this and replace with trial 51. trial looks fine in Nexus though...
    end
    
    % Number of trials in the dataset:
    num_trials = length(TrialData);
    distA = nan*ones(size(num_trials)); distL = nan*ones(size(num_trials));
    
    % Get mass of subj for ang momentum metric. Since compare solo to
    % assist beam, only care subj's with dynamics data
    if subj == 3
        mass = 49.3; % kg
        height = 1.60; % m
    elseif subj == 4
        mass = 56.9; % kg
        height = 1.61; % m
    elseif subj == 5
        mass = 73.8; % kg
        height = 1.64; % m
    elseif subj == 6
        mass = 78.2; % kg
        height = 1.73; % m
    elseif subj == 7
        mass = 76.0; % kg
        height = 1.82; % m
    elseif subj == 8
        mass = 67.5; % kg
        height = 1.69; % m
    elseif subj == 9
        mass = 75.2; % kg
        height = 1.68; % m
    elseif subj == 10
        mass = 55.8; % kg
        height = 1.61; % m
    elseif subj == 11
        mass = 62.0; % kg
        height = 1.64; % m
    elseif subj == 12
        mass = 67.2; % kg
        height = 1.72; % m
    elseif subj == 13
        mass = 54.8; % kg
        height = 1.61; % m
    elseif subj == 14
        mass = 69.3; % kg
        height = 1.78; % m
    end
        
    clear Markers
    temp = [];
    firstPlot = 0; % initialize figures if first time in plotCheck

    %% Trial Analysis Main Loop
    for n = 1:num_trials
        clear start_idx stop_idx APMarkers LHEE RHEE RTOE LPSI RPSI beamMidline
        disp(['Processing ',TrialData(n).Info.Trial]);
        trial = str2num(TrialData(n).Info.Trial(end-1:end));
        sample_rate = TrialData(n).Markers.samplerate;
        
        % Check for trials to process
        if any(strcmpi(TrialData(n).Info.Condition, to_process))
            %% Get markers from each trial type %%
            % Median filter all of the marker data to remove jumps. Then
            % calculate the start and stop indices for HHI analysis
            max_distance = TrialData(n).Info.Distance_Traveled.*25.4; %convert from in to mm
            if strcmpi(TrialData(n).Info.Condition, 'Assist Solo') % Assist Solo
                Markers = TrialData(n).Markers;
                beamTrial = 0;
                if isfield(Markers,'AP') == 1 % Some trials denoted assistant by AP while others did not
                    Markers.AP = medfiltFields(Markers.AP,1);
                    TimeMarkers = Markers.AP;
                    APMarkers = Markers.AP;
                    LHEE = Markers.AP.LHEE;
                    RHEE = Markers.AP.RHEE;
                    RTOE = Markers.AP.RTOE;
                    LPSI = Markers.AP.LPSI;
                    RPSI = Markers.AP.RPSI;
                else % Trials not denote assistant as "AP"
                    Markers = medfiltFields(TrialData(n).Markers,1);
                    TimeMarkers = Markers;
                    APMarkers = Markers;
                    LHEE = Markers.LHEE;
                    RHEE = Markers.RHEE;
                    RTOE = Markers.RTOE;
                    LPSI = Markers.LPSI;
                    RPSI = Markers.RPSI;
                end
            elseif any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist Ground','Assist beam'}))
                % Assist Ground, Assist Beam
                Markers = TrialData(n).Markers;
                Markers.AP = medfiltFields(Markers.AP,1);
                Markers.POB = medfiltFields(Markers.POB,1);
                TimeMarkers = Markers.POB;
                APMarkers = Markers.AP;
                LHEE = Markers.POB.LHEE;
                RHEE = Markers.POB.RHEE;
                RTOE = Markers.POB.RTOE;
                LPSI = Markers.POB.LPSI;
                RPSI = Markers.POB.RPSI;
            else % Solo
                % Solo Beam, Solo Ground
                Markers = medfiltFields(TrialData(n).Markers,1);
                TimeMarkers = Markers;
                LHEE = Markers.LHEE;
                RHEE = Markers.RHEE;
                RTOE = Markers.RTOE;
                LPSI = Markers.LPSI;
                RPSI = Markers.RPSI;
            end
            
            %% Time window of analysis
            
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist beam','Solo Beam','Solo beam'}))
                beamTrial = 1;
            else
                beamTrial = 0;
            end
            
            % Torso pos used for finding start/stop of trial
            clear torso torsoFilt vTorso vTorsoFilt aTorso aTorsoFilt
            % Get torso values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
            if (subj == 4 && trial == 7) || (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                torso = TimeMarkers.C7./1000;
            else
                torso = TimeMarkers.CLAV./1000;
            end
            
            vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
            
            % Luke's method using max_distance
            [start_L,stop_L] = getHHIAnalysisWindow(TimeMarkers,max_distance);
            % Check the assister's arm
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist beam', 'Assist Ground', 'Assist Solo'}))
                [start_L,stop_L] = checkArm(APMarkers,start_L,stop_L);
            end

            % New method
            [start_idx,stop_idx,temp.ht,temp.midline] = getHHIAnalysisWindow_MW(TimeMarkers,vLHEEZfilt,beamTrial,torso,subj,trial);
            TrialData(n).Results.ht = temp.ht; % (m)
            TrialData(n).Results.midline = temp.midline; % (m)
            
            % Store torso state for later analysis/plotting
            TrialData(n).Results.torso = torso(start_idx:stop_idx,:); % (m)
            for k = 1:3
                torsoFilt(:,k) = filtfilt(sos, g, torso(:,k)); % (m) % Use only for taking deriv.
            end
            vTorso = diff(torsoFilt).*sample_rate; 
            % Filter vel before take derivative. 
            for k = 1:3
                vTorsoFilt(:,k) = filtfilt(sos, g, vTorso(:,k));
            end
            aTorso = diff(vTorsoFilt).*sample_rate;
            % Filter acc to use for regression later
            for k = 1:3
                aTorsoFilt(:,k) = filtfilt(sos, g, aTorso(:,k));
            end
            % Trim data to that which corresponds to force data
            vTorsoFilt = vTorsoFilt(start_idx+1:stop_idx+1,:);
            aTorsoFilt = aTorsoFilt(start_idx+2:stop_idx+2,:);
            
            % Store torso derivatives for later analysis/plotting
            TrialData(n).Results.vTorso = vTorso(start_idx+1:stop_idx+1,:);
            TrialData(n).Results.aTorso = aTorso(start_idx+2:stop_idx+2,:);
            
            %% Construct the time axis
            clear time
            time = 0:(stop_idx - start_idx);
            time = time./sample_rate;
            
            %% Get POB COM position as avg L and R PSIS % don't use this for now bc lots of work to check marker fill correct for both PSIS markers all trials during start-stop period
            COM = (LPSI + RPSI)/(2*1000); % Take mean
            
            %% Calculate kinematic metrics 
            % Use the medfiltered data from Markers struct
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam','Assist beam','Solo Beam','Solo beam','Assist Ground', 'Solo Ground'}))
                %% Lateral sway and hip and leg angles
                sway = torso(start_idx:stop_idx,1);
                if any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam', 'Solo Ground'}))                   
                    % Also collect thorax (approx with
                    % projection of line between two thorax markers) ang
                    % rotation for POB. Calculate leg segment obliq by
                    % 1) project all malleolus markers and PSIS markers
                    % into frontal plane, 2) take midpoint 2 malleoli, 3)
                    % take midpoint PSIS markers, 4) take vector between
                    % the two points and get angle in frontal plane
                    TrialData(n).Results.pelvicObliq = TrialData(n).Markers.LPelvisAngles(start_idx:stop_idx,1)./1000; % L and R pelvis angles are exact mirror images of each other. Pelvic data not really necessary, just look at thorax and leg segments for now     
                    % Leg segment angle (only care x and z components since
                    % will look at frontal plane angle
                    temp.mANK(:,1) = mean([Markers.LANK(start_idx:stop_idx,1) Markers.RANK(start_idx:stop_idx,1)],2)./1000;
                    temp.mANK(:,3) = mean([Markers.LANK(start_idx:stop_idx,3) Markers.RANK(start_idx:stop_idx,3)],2)./1000;
                    temp.mPSIS(:,1) = mean([Markers.LPSI(start_idx:stop_idx,1) Markers.RPSI(start_idx:stop_idx,1)],2)./1000;
                    temp.mPSIS(:,3) = mean([Markers.LPSI(start_idx:stop_idx,3) Markers.RPSI(start_idx:stop_idx,3)],2)./1000;
                    temp.legVec = temp.mPSIS-temp.mANK;
                    temp = [];
                else % Force data exists, two people marker set
                    % Also collect thorax (approx with
                    % projection of line between two thorax markers) ang
                    % rotation for POB
                    if isfield(Markers.POB,'LPelvisAngles')
                        TrialData(n).Results.pelvicObliq = TrialData(n).Markers.POB.LPelvisAngles(start_idx:stop_idx,1); % L and R pelvis angles are exact mirror images of each other  
                    else
                        TrialData(n).Results.pelvicObliq = nan;
                    end
                    % Leg segment angle (only care x and z components since
                    % will look at frontal plane angle
                    temp.mANK(:,1) = mean([Markers.POB.LANK(start_idx:stop_idx,1) Markers.POB.RANK(start_idx:stop_idx,1)],2)./1000;
                    temp.mANK(:,3) = mean([Markers.POB.LANK(start_idx:stop_idx,3) Markers.POB.RANK(start_idx:stop_idx,3)],2)./1000;
                    temp.mPSIS(:,1) = mean([Markers.POB.LPSI(start_idx:stop_idx,1) Markers.POB.RPSI(start_idx:stop_idx,1)],2)./1000;
                    temp.mPSIS(:,3) = mean([Markers.POB.LPSI(start_idx:stop_idx,3) Markers.POB.RPSI(start_idx:stop_idx,3)],2)./1000;
                    temp.legVec = temp.mPSIS-temp.mANK;
                    temp = [];
                end
                temp = [];
                
                sway = sway - TrialData(n).Results.midline; % same as torso pos with midline removed
                TrialData(n).Results.beamerSway = sway;
                TrialData(n).Results.beamerCOMSway = COM(start_idx:stop_idx,1);
                
                %% Distance is calculated as the difference between the
                % forward positions of the Rheel at the beginning of the
                % analysis window and the heel of the forward foot at the
                % end of the analysis window. Do only for beam trials!
                if any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam', 'Assist Beam','Assist beam'}))
                    dist = torso(end,2) - torso(1,2); % torso is already cut to use start and stop indices
                    if dist >= 3.65 % length of beam
                        distA(n) = 3.65;
                    else
                        distA(n) = dist; 
                    end
                    TrialData(n).Results.totalDistance = distA(n); 
                    %% Luke's method on finding beam distance
                    if any(strcmpi(TrialData(n).Info.Condition,'Solo Beam'))
                        torsoY = TrialData(n).Markers.CLAV(:,2);
                    else
                        torsoY = TrialData(n).Markers.POB.CLAV(:,2);
                    end
                    tempD = torsoY(stop_L) - torsoY(start_L);
                    if (TrialData(n).Info.Distance_Traveled == 144)
                        distL(n) = 144*25.4;
                    else
                        distL(n) = tempD;
                    end
                else
                    TrialData(n).Results.totalDistance = nan;
                    distA(n) = nan; distL(n) = nan;
                end
                
                %% Angular state and momentum of POB
                clear angTorso wTorso alphaTorso I Ly
                angTorso = atan2(torso(:,1)-TrialData(n).Results.midline,torso(:,3)); % positive CW, relative to vertical
                wTorso = diff(filtfilt(sos,g,angTorso)).*sample_rate; 
                alphaTorso = diff(filtfilt(sos,g,wTorso)).*sample_rate;
                TrialData(n).Results.angTorso = angTorso(start_idx:stop_idx);
                TrialData(n).Results.wTorso = wTorso(start_idx+1:stop_idx+1);
                TrialData(n).Results.alphaTorso = alphaTorso(start_idx+2:stop_idx+2);
                I = mass*height^2/3;
                TrialData(n).Results.Ly = I*TrialData(n).Results.wTorso; 
                
                %% We also store the average speed of the trial, and a
                % logical value indicating if the trial was completed.
                TrialData(n).Results.avgSpeed = TrialData(n).Results.totalDistance./time(end); % time has been trimmed to cover start to stop only
                TrialData(n).Results.completed = (TrialData(n).Info.Distance_Traveled == 144);

            elseif strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                TrialData(n).Results.vTorso = nan;
                TrialData(n).Results.aTorso = nan;
                TrialData(n).Results.angTorso = nan;
                TrialData(n).Results.wTorso = nan;
                TrialData(n).Results.alphaTorso = nan;
                TrialData(n).Results.Ly = nan;
                TrialData(n).Results.totalDistance = nan;
                TrialData(n).Results.avgSpeed = nan;
                TrialData(n).Results.pelvicObliq = nan;
                distA(n) = nan; distL(n) = nan;
                TrialData(n).Results.beamerSway = nan;
            end
            
            %% Calculations depending on force data, 2 person marker set
            if ~isempty(strfind(TrialData(n).Info.Condition, 'Assist')) % Check there is an assistant
                % Do some of this analysis for Assist Solo
                %% Window and Transpose Force and Marker Data %% 
                % Now that we have the start and stop indices, we can get the
                % force data and the marker data for the time when the subject
                % is walking.               
                
                % Get the sampling rate for the forces                
                force_rate = TrialData(n).Forces.samplerate;
                
                clear Force Torque F T
                % Collect the force and torque voltages into the initialized arrays
                for k = 1:3 % for each direction                  
                    % Downsample the forces and torques to be at the same
                    % sampling rate as the marker data
                    F = resample(TrialData(n).Forces.(force_names{k}), sample_rate, force_rate);
                    T = resample(TrialData(n).Forces.(torque_names{k}), sample_rate, force_rate);
                    % Median filter forces to reduce noise - why do this
                    % after downsampling? Perhaps to be similar to
                    % processing of marker data?
                    F = medfilt1(F);
                    T = medfilt1(T);
                    % Limit the forces and torques to the same analysis window
                    % as the marker data
                    Force(:,k) = F(start_idx:stop_idx);
                    Torque(:,k) = T(start_idx:stop_idx); % This is torque on sensor or IP not on POB COM
                end
                
                if isfield(Markers,'FH')
                    % Median filter force handle MARKERS to remove jumps
                    FH = medfiltFields(TrialData(n).Markers.FH,1); % Tolerance check in later recoverForces code looks at mm so don't convert to m
                    % Get only the markers in the trial window
                    FH = windowMarkers(FH,start_idx,stop_idx);
                    % Transpose the force handle markers for vector
                    % compatibility
                    FH = transposeMarkers(FH);
                    % Save raw forces to check biases (this is after
                    % changing signs on Fx and Fy voltage
                    TrialData(n).Results.FV = Force;
                    % Recover forces in the VICON frame
                    Force = recoverForces(Force',Torque',FH); 
                    Force = Force';
                end
                
                % Up to now all forces are in sensor's coordinate frame,
                % change it to the POB's coordinate frame so it's force
                % felt by person. This will then be consistent for the
                % model that fits torque felt by person (on the person) to
                % torso state. Do this after 
                % recoverForces projects force signals into lab CS since
                % that code subtracts off a signed baseline voltage
                % For sensor, tension > 0 for all voltages. 
                Force = -Force;
                % Store the data
                TrialData(n).Results.Force = Force;
                
                % Calculate 2D vector norm of force
                clear Fmag theta
                for i = 1:length(Force(:,1))
                    Fmag(i) = norm(Force(i,[1 3]));
                    theta(i) = -(atan2(Force(i,3),Force(i,1))-pi/2); % CW > 0, relative to vertical
                end
                               
                % Need to subtract out force of gravity on FH for vertical
                % dir and subtract out force of accelerating mass of FH in
                % each direction (?)              
                
                %% --------------- Dynamics ANALYSIS --------------------%%
                if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo') % Do only for Assist
                    %% Get markers used for power calculations
                    if isfield(Markers,'AP') == 1                   
                        RSHO = Markers.AP.RSHO(start_idx:stop_idx, :)./1000;
                        RFIN = Markers.AP.RFIN(start_idx:stop_idx, :)./1000;
                    else
                        RSHO = Markers.RSHO(start_idx:stop_idx, :)./1000;
                        RFIN = Markers.RFIN(start_idx:stop_idx, :)./1000;
                    end

                    if isfield(Markers,'POB') == 1
                        LSHO = Markers.POB.LSHO(start_idx:stop_idx, :)./1000;
                        LFIN = Markers.POB.LFIN./1000; 
                    else
                        LSHO = Markers.LSHO(start_idx:stop_idx, :)./1000;
                        LFIN = Markers.LFIN./1000;
                    end

                    % Filter kinematics and force
                    clear ForceFilt rFin lFin rFilt lFilt vFin vFinFilt RASIFilty armAsst armPOB POBpower IP_power_XZ
                    for k = 1:3
                        ForceFilt(:,k) = filtfilt(sos, g, Force(:,k));
                        rFilt(:,k) = filtfilt(sos, g, RFIN(:,k)); 
                        lFilt(:,k) = filtfilt(sos, g, LFIN(:,k)); 
                        vFin(:,k) = diff(lFilt(:,k)).*sample_rate;
                    end
                    rFin = RFIN; % Don't double-filter! Use filtered pos only for taking derivatives and doing regression and corr later
                    lFin = LFIN(start_idx:stop_idx,:);
                    
%                   vFin = diff(rFilt).*sample_rate;
                    TrialData(n).Results.IP = lFin; 
                    for i = 1:3
                        vFinFilt(:,i) = filtfilt(sos, g, vFin(:,i));
                        aFin = diff(vFinFilt).*sample_rate; 
                    end
 
                    % Trim data to that which corresponds to force
                    vFin = vFin(start_idx+1:stop_idx+1,:);
                    vFinFilt = vFinFilt(start_idx+1:stop_idx+1,:);
                    aFin = aFin(start_idx+1:stop_idx+1,:);
                    
                    %% Power exchange at IP in x, z, and xz vector dir's
                    clear IP_power IP_powerX_hi IP_powerX_lo IP_cumWork  IP_netWork_norm IP_power_XZ 
                    % Don't filter for power calc's
                    for i = 1:3
                        IP_power(:,i) = Force(:,i).*vFin(:,i); % P < 0 dissipates energy (F < 0 for tension on POB and v > 0 for move to the right)
                        if i == 1 % x axis, check worst case drift effect of 1.5N
                            IP_powerX_hi(1,:) = (Force(:,i)+1.5).*vFin(:,i);
                            IP_powerX_lo(1,:) = (Force(:,i)-1.5).*vFin(:,i);
                        end
                        IP_cumWork(:,i) = cumsum(IP_power(:,i))./sample_rate;
                        IP_netWork_norm(i) = sum(IP_power(:,i))/sample_rate/(time(end)-time(1)); % One number characterizing whole trial overall, normalize by time length of trial
                    end
                    
                    % Power in 2D vector's direction
                    for j = 1:length(vFin)
                        IP_power_XZ(j) = dot(Force(j,[1 3]),vFin(j,[1 3]));
                    end
                    
                    %% Power on POB in x, z, and xz vector dir's
                    clear POB_power POB_powerX_hi POB_powerX_lo POB_cumWork  POB_netWork_norm POB_power_XZ 
                    % Don't filter for power calc's
                    for i = 1:3
                        POB_power(:,i) = Force(:,i).*TrialData(n).Results.vTorso(:,i); % P < 0 dissipates energy (F < 0 for tension on POB and v > 0 for move to the right)
                        if i == 1 % x axis, check worst case drift effect of 1.5N
                            POB_powerX_hi(1,:) = (Force(:,i)+1.5).*TrialData(n).Results.vTorso(:,i);
                            POB_powerX_lo(1,:) = (Force(:,i)-1.5).*TrialData(n).Results.vTorso(:,i);
                        end
                        POB_cumWork(:,i) = cumsum(POB_power(:,i))./sample_rate;
                        POB_netWork_norm(i) = sum(POB_power(:,i))/sample_rate/(time(end)-time(1)); % One number characterizing whole trial overall, normalize by time length of trial
                    end
                    
                    % Power in 2D vector's direction
                    for j = 1:length(TrialData(n).Results.vTorso)
                        POB_power_XZ(j) = dot(Force(j,[1 3]),TrialData(n).Results.vTorso(j,[1 3]));
                    end
                    
                    %% Calculate torque and ang power on POB. Take xz vector for moment arm and force
                    % Examine torque about torso, COM, and beam (heel)
                    clear rTorso rCOM theta rTorsoxz rGdxz rCOMht
                    rTorso = lFin - torso(start_idx:stop_idx,:);
%                     rTorso = rFin - torso(start_idx:stop_idx,:);
    %                 rCOM = rFin - COM(start_idx:stop_idx,:); % Need to check
    %                 COM markers ID'ed correctly before do this calculation
                    if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist beam'})) % Use beam midline
                        rGd = lFin - [TrialData(n).Results.midline*ones(length(lFin),1) torso(start_idx:stop_idx,2) TrialData(n).Results.ht*ones(length(lFin),1)]; 
%                         rGd = rFin - [TrialData(n).Results.midline*ones(length(rFin),1) torso(start_idx:stop_idx,2) TrialData(n).Results.ht*ones(length(rFin),1)]; 
                    else % Overgound, assume axis is on ground at z = 0
                        rGd = lFin - [TrialData(n).Results.midline*ones(length(lFin),1) torso(start_idx:stop_idx,2) zeros(length(lFin),1)]; % Use midline calculated from mean torso lateral pos during period of interest
%                         rGd = rFin - [TrialData(n).Results.midline*ones(length(rFin),1) torso(start_idx:stop_idx,2) zeros(length(rFin),1)]; % Use midline calculated from mean torso lateral pos during period of interest
                    end

                    % Find angle between two vec's in xz plane to keep for
                    % reference. Also keep moment arm         
                    for i = 1:length(rGd)
                        rTorsoxz(i) = norm(rTorso(i,[1 3]));
                        rGdxz(i) = norm(rGd(i,[1 3]));
                    end
                    theta = -atan2(Force(:,3),Force(:,1)) + atan2(rTorso(:,3),rTorso(:,1));
                    TrialData(n).Results.maTorso = rTorsoxz'.*sin(theta);
                    theta = -atan2(Force(:,3),Force(:,1)) + atan2(rGd(:,3),rGd(:,1));
                    TrialData(n).Results.maGd = rGdxz'.*sin(theta);

                    % Get torque from 3D cross-product and take My component
                    % Also calc worst case drift scenario
                    temp = cross(rTorso,Force);
                    TrialData(n).Results.TyTorso = temp(:,2); % Torque about y axis
    %                 TrialData(n).Results.torqueCOM = cross(rCOM,Force);               
                    temp = cross(rGd,Force);               
                    TrialData(n).Results.TyGd = temp(:,2); % Net torque on POB
                    clear r TyFx TyFz temp
                    TrialData(n).Results.TyFx = rGd(:,3).*TrialData(n).Results.Force(:,1);
                    TrialData(n).Results.TyFz = -rGd(:,1).*TrialData(n).Results.Force(:,3); % Must change sign here to make torque sign correct
                    
                    TrialData(n).Results.TyInt = calcTint(TrialData(n).Results.TyFx,TrialData(n).Results.TyFz);
                    % Compute torque due to gravity acting at 1/2 of height
                    % of POB
                    rCOMht = height/2;
                    TrialData(n).Results.Tg = rCOMht*sin(TrialData(n).Results.angTorso)*9.81*mass;
                    
                    TrialData(n).Results.TyExt = TrialData(n).Results.TyGd - TrialData(n).Results.Tg; % Is external torque (applied by partner) the difference between net torque and gravitational torque?
%                     % Check sign convention
%                     ind = find(temp(:,2) > 0,1,'first');
%                     quiver(0,0,rGd(ind,1),rGd(ind,3)),hold on,quiver(rGd(ind,1),rGd(ind,3),Force(ind,1),Force(ind,3));
                    
                    % Effects of worst case drift on TyFx, TyFz, and TyInt 
                    
                    TrialData(n).Results.TyFxh = rGd(:,3).*(TrialData(n).Results.Force(:,1)+1.5);
                    TrialData(n).Results.TyFxl = rGd(:,3).*(TrialData(n).Results.Force(:,1)-1.5);
                    TrialData(n).Results.TyFzh = -rGd(:,1).*(TrialData(n).Results.Force(:,3)+5.7);
                    TrialData(n).Results.TyFzl = -rGd(:,1).*(TrialData(n).Results.Force(:,3)-5.7);
                    TrialData(n).Results.TyIntXhZh = calcTint(TrialData(n).Results.TyFxh,TrialData(n).Results.TyFzh);
                    TrialData(n).Results.TyIntXlZh = calcTint(TrialData(n).Results.TyFxl,TrialData(n).Results.TyFzh);
                    TrialData(n).Results.TyIntXhZl = calcTint(TrialData(n).Results.TyFxh,TrialData(n).Results.TyFzl);
                    TrialData(n).Results.TyIntXlZl = calcTint(TrialData(n).Results.TyFxl,TrialData(n).Results.TyFzl);
                    % Power
                    TrialData(n).Results.angP = (TrialData(n).Results.TyGd).*TrialData(n).Results.wTorso;
                    TrialData(n).Results.angPFx = (TrialData(n).Results.TyFx).*TrialData(n).Results.wTorso;
                    TrialData(n).Results.angPFz = (TrialData(n).Results.TyFz).*TrialData(n).Results.wTorso;

                    %% Calculate xcorr for a variety of signals to compare to LT work and to adjust regression if needed
                    % for xcorr(x,y), y lags x if lag < 0;

                    % Use filtered signals before xcorr so corr is not driven by noise
                    
                    % ML/x dir
                    % F_IP to POB torso disp
                    temp.swayX = torsoFilt(start_idx:stop_idx,1) - nanmean(torsoFilt(start_idx:stop_idx,1)); % This is not displacement! This is position with mean removed!
                    [temp.a,temp.b] = getXcorr(ForceFilt(:,1),temp.swayX,sample_rate);
                    TrialData(n).Results.lagFIPTorsoX = temp.a;
                    TrialData(n).Results.xcorrFIPTorsoX = temp.b;
                    temp = [];
                  
                    % F_IP to POB vTorso in ML dir only       
                    [temp.a,temp.b] = getXcorr(ForceFilt(:,1),vTorsoFilt(:,1),sample_rate);
                    TrialData(n).Results.lagFIPvTorsoX = temp.a;
                    TrialData(n).Results.xcorrFIPvTorsoX = temp.b;
                    temp = [];
                    
                    % vIP to vTorso in ML dir only. 
                    [temp.a,temp.b] = getXcorr(vFinFilt(:,1),vTorsoFilt(:,1),sample_rate);
                    TrialData(n).Results.lagvIPvTorsoX = temp.a;
                    TrialData(n).Results.xcorrvIPvTorsoX = temp.b;
                    temp = [];
                    
                    % Vert/z dir
                    % F_IP to POB torso disp
                    temp.pos = torsoFilt(start_idx:stop_idx,3) - nanmean(torsoFilt(start_idx:stop_idx,3)); % This is not displacement! This is position with mean removed!
                    [temp.a,temp.b] = getXcorr(ForceFilt(:,3),temp.pos,sample_rate);
                    TrialData(n).Results.lagFIPTorsoZ = temp.a;
                    TrialData(n).Results.xcorrFIPTorsoZ = temp.b;
                    temp = [];
                  
                    % F_IP to POB vTorso in ML dir only       
                    [temp.a,temp.b] = getXcorr(ForceFilt(:,3),vTorsoFilt(:,3),sample_rate);
                    TrialData(n).Results.lagFIPvTorsoZ = temp.a;
                    TrialData(n).Results.xcorrFIPvTorsoZ = temp.b;
                    temp = [];
                    
                    % vIP to vTorso in ML dir only. 
                    [temp.a,temp.b] = getXcorr(vFinFilt(:,3),vTorsoFilt(:,3),sample_rate);
                    TrialData(n).Results.lagvIPvTorsoZ = temp.a;
                    TrialData(n).Results.xcorrvIPvTorsoZ = temp.b;
                    temp = [];
                    
                    % Ty corr to torso state
                    % ang disp
                    T = filtfilt(sos, g, TrialData(n).Results.TyGd);
                    temp.ang = filtfilt(sos, g, TrialData(n).Results.angTorso);
                    [temp.a,temp.b] = getXcorr(T,temp.ang,sample_rate);
                    TrialData(n).Results.lagTyAngTorso = temp.a;
                    TrialData(n).Results.xcorrTyAngTorso = temp.b;
                    temp = [];

                    % ang vel
                    temp.w = filtfilt(sos,g,TrialData(n).Results.wTorso);
                    [temp.a,temp.b] = getXcorr(T,temp.w,sample_rate);
                    TrialData(n).Results.lagTyWTorso = temp.a;
                    TrialData(n).Results.xcorrTyWTorso = temp.b;
                    temp = [];    
                    
                    clear T
                    
                    %% Fit F_IP to torso state var's 

                    % Regress per direction. Regress
                    % iteratively, throwing away any variables with coeff's
                    % whose CI's include zero. For ML disp, look at torso
                    % relative to midline. Assume model acts to
                    % restore torso to midline in ML dir and to torso
                    % height and sagittal pos at beg of trial. Must reverse
                    % sign of force if want positive damper and spring coeff
                    % to mean passive elements. 
                    % Up to this point, sign convention for all force
                    % directions is positive when tension or shear for POB
                    % away from IP. For fit model, regress to -m is consistent with a
                    % positive spring or damper restoring person to midline
                    % or starting position
                    
                    % X/ML dir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clear m;
                    lag = floor(.157*sample_rate); % calculated from xcorr analysis
                    i = 1; % Just do for x dir. Not make sense to fit model to AP/y dir since no equilibrium position in y dir.
                    ref = TrialData(n).Results.midline; 
                    % No lag
                    m = [TrialData(n).Results.aTorso(:,i) TrialData(n).Results.vTorso(:,i) TrialData(n).Results.torso(:,i)-ref ones(size(ForceFilt(:,i)))]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                    [cx_torso,rsqx_torso,px_torso,cx_torso_st,rsqx_torso_st] = regressIterNew(ForceFilt(:,i),-m); % Fx > 0 if tension. If tension and displacement > 0 (POB move to R), should get positive spring
                                        
                    % Lag of 157ms (kinem's lag force?)  - check this code
%                     mlag = [aTorso(start_idx-2+lag:stop_idx-2+lag,i) vTorso(start_idx-1+lag:stop_idx-1+lag,i) torso(start_idx+lag+2:end,i)-ref ones(size(Force(lag:end,i)))]; 
%                     [cx_lag_torso,rsqx_lag_torso,px_lag_torso] = regressIter(Force(lag+2:end,i),mlag); 
%                     [cx_lag_torso,rsqx_lag_torso,px_lag_torso,cx_lag_torso_st,rsqx_lag_torso_st] = regressIterNew(Force(:,i),mlag);
                    % Test correlations
                    % no lag
%                     TrialData(n).Results.FresX = Force(:,i) - m*cx_torso;  % Calculate corr power to residual. However, oppositte force sign conventions for power vs. force for model, so just know that -rho values mean correlated.
%                     [temp, TrialData(n).Results.corrPFxpval] = corr(TrialData(n).Results.FresX,IP_power(:,i));
%                     if TrialData(n).Results.corrPFxpval < 0.05
%                         TrialData(n).Results.corrPowerFresX = temp;
%                     else
%                         TrialData(n).Results.corrPowerFresX = nan;
%                     end
                    % lag
%                     TrialData(n).Results.FresLagX = Force(lag+2:end,i) - mlag*cx_lag_torso;  % Calculate corr power to residual. However, oppositte force sign conventions for power vs. force for model, so just know that -rho values mean correlated.
%                     [temp, TrialData(n).Results.corrPFxLagpval] = corr(TrialData(n).Results.FresLagX,IP_power(:,i));
%                     if TrialData(n).Results.corrPFxLagpval < 0.05
%                         TrialData(n).Results.corrPowerFresLagX = temp;
%                     else
%                         TrialData(n).Results.corrPowerFresLagX = nan;
%                     end           
                    
                    % Z/Vert dir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clear m;
                    i = 3; 
                    eqbmHt = TrialData(n).Results.torso(1,i); % Use height of torso right after step onto beam as position reference for fit model in vertical/z dir. What we want is height of torso when person is upright and on balance...
%                     eqbmHt = rFin(1,3); 
                    TrialData(n).Results.eqbmHt = eqbmHt; % Same for every trial
                    % No lag
                    m = [TrialData(n).Results.aTorso(:,i) TrialData(n).Results.vTorso(:,i) TrialData(n).Results.torso(:,i)-eqbmHt ones(size(ForceFilt(:,i)))]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                    [cz_torso,rsqz_torso,pz_torso,cz_torso_st,rsqz_torso_st] = regressIterNew(ForceFilt(:,i),-m); % Fz > 0 if tension. If tension and displacement > 0 , should get positive spring
                     
                    % 2D xz vector dir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clear m;
                    i = [1 3]; 
                    eqbm = [ref eqbmHt]; % Use ref's from X and Z regressions above
                    % No lag
                    % Use dot product to find scalar projection of torso state
                    % vectors onto force vector
                    for k = 1:length(ForceFilt(:,1))
                        temp.F = ForceFilt(k,i);
                        temp.Fn(k) = norm(temp.F);
                        temp.a(k) = dot(TrialData(n).Results.aTorso(k,i),temp.F)/temp.Fn(k);
                        temp.v(k) = dot(TrialData(n).Results.vTorso(k,i),temp.F)/temp.Fn(k);
                        temp.d(k) = dot((TrialData(n).Results.torso(k,i)-eqbm),temp.F)/temp.Fn(k);
                    end
                    m = [temp.a' temp.v' temp.d' ones(size(temp.a'))]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                    [cxz_torso,rsqxz_torso,pxz_torso,cxz_torso_st,rsqxz_torso_st] = regressIterNew(temp.Fn',-m); % Fz > 0 if tension. If tension and displacement > 0 , should get positive spring
                    
%                     % Plot check projection of each vector is correct at a
%                     % given instance in time
%                     k = 100;
%                     quiver(0,0,ForceFilt(k,1),ForceFilt(k,3)),hold on;
%                     quiver(0,0,TrialData(n).Results.aTorso(k,1),TrialData(n).Results.aTorso(k,3));
%                     quiver(0,0,temp.a(k)*ForceFilt(k,1)/temp.Fn(k),temp.a(k)*ForceFilt(k,3)/temp.Fn(k));
%                     legend('F','a','a proj');
%                     
%                     quiver(0,0,ForceFilt(k,1),ForceFilt(k,3)),hold on;
%                     quiver(0,0,TrialData(n).Results.vTorso(k,1),TrialData(n).Results.vTorso(k,3));
%                     quiver(0,0,temp.v(k)*ForceFilt(k,1)/temp.Fn(k),temp.v(k)*ForceFilt(k,3)/temp.Fn(k));
%                     legend('F','v','v proj');
%                     
%                     quiver(0,0,ForceFilt(k,1),ForceFilt(k,3)),hold on;
%                     quiver(0,0,TrialData(n).Results.torso(k,1)-eqbm(1),TrialData(n).Results.torso(k,3)-eqbm(2));
%                     quiver(0,0,temp.d(k)*ForceFilt(k,1)/temp.Fn(k),temp.d(k)*ForceFilt(k,3)/temp.Fn(k));
%                     legend('F','disp','disp proj');
                    
                    temp = [];
                    
                    %% Fit Ty Gd to torso ang state var's. Use filtered versions of all var's
                    
                    temp.Tnet = filtfilt(sos,g,TrialData(n).Results.TyGd);
                    temp.TFx = filtfilt(sos,g,TrialData(n).Results.TyFx);
                    temp.TFz = filtfilt(sos,g,TrialData(n).Results.TyFz);
                    temp.Tint = filtfilt(sos,g,TrialData(n).Results.TyInt);

                    % Regress iteratively, throwing away any variables with coeff's
                    % whose CI's include zero. For torsional spring, look at torso
                    % relative to angle of pi/2 (vertical). Want sign
                    % convention of positive coefficient means impeding or
                    % passive element so need to regress to -m
                                       
                    clear m;
                    ref = 0; 
                    % No lag
                    temp.theta = filtfilt(sos,g,angTorso);
                    angTorsoFilt = temp.theta(start_idx:stop_idx);
                    temp.w = filtfilt(sos,g,wTorso);
                    wTorsoFilt = temp.w(start_idx+1:stop_idx+1);
                    temp.alpha = filtfilt(sos,g,alphaTorso);
                    alphaTorsoFilt = temp.alpha(start_idx+2:stop_idx+2);
                                        
                    m = [alphaTorsoFilt wTorsoFilt angTorsoFilt-ref ones(size(temp.Tnet))]; 
                    [c_ang_torso,rsq_ang_torso,p_ang_torso,c_ang_torso_st,rsq_ang_torso_st] = regressIterNew(temp.Tnet,-m);                     
                    [c_angFx_torso,rsq_angFx_torso,p_angFx_torso,c_angFx_torso_st,rsq_angFx_torso_st] = regressIterNew(temp.TFx,-m);                     
                    [c_angFz_torso,rsq_angFz_torso,p_angFz_torso,c_angFz_torso_st,rsq_angFz_torso_st] = regressIterNew(temp.TFz,-m);    
                    
                    clear m;                     
                    m = [abs(alphaTorsoFilt) abs(wTorsoFilt) abs(angTorsoFilt-ref) ones(size(temp.Tint))]; % Use abs values of state b/c logic is that mag of interact torque increases as deviation from vertical position and zero speed/acc increases
                    [c_angInt_torso,rsq_angInt_torso,p_angInt_torso,c_angInt_torso_st,rsq_angInt_torso_st] = regressIterNew(temp.Tint,-m);                        
                    
                    %% Get POB RASI marker to calculate AP velocity of whole body 
                    RASIFilty = filtfilt(sos, g, Markers.POB.RASI(start_idx:stop_idx,2))./1000; % (m)
                    % Store pos for later analysis/plotting
                    TrialData(n).Results.RASI = Markers.POB.RASI(start_idx:stop_idx,2)./1000;

                    % Need to filter first before taking derivative!
                    % Store velocity for later analysis/plotting
                    TrialData(n).Results.vyRASI = diff(RASIFilty).*sample_rate;
                    
                    %% Get Assistant's arm len state (vector) and POB's arm len state (vector)
                    % Compute the Arm Lengths as just scalar magnitude 
                    for j = 1:length(RFIN)
                        armAsst(j) = norm(RFIN(j,:) - RSHO(j,:)); % (m)
                        if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                            armPOB(j) = norm(LSHO(j,:) - LFIN(j,:)); % (m)
                        end
                    end
                    % Filter before take derivative
                    armAsst = filtfilt(sos, g, armAsst);
                    if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                        armPOB = filtfilt(sos, g, armPOB);
                    end
                    % vel
                    vArmAsst = diff(armAsst).*sample_rate;
                    if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                        vArmPOB = diff(armPOB).*sample_rate;
                    end
                    % Filter vel before take derivative
                    vArmAsstFilt = filtfilt(sos, g, vArmAsst);
                    if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                        vArmPOBFilt = filtfilt(sos, g, vArmPOB);
                        aArmPOB = diff(vArmPOBFilt).*sample_rate;
                    end
                    aArmAsst = diff(vArmAsstFilt).*sample_rate;
                    
                    %% Get POB's effective arm length only in ML dir
                    TrialData(n).Results.armPOBX = RFIN(:,1) - RSHO(:,1);
                    
                end
                
%% --------------- STORE THE DYNAMICS  RESULTS ----------------------------------- %%
                % For assisted trials
                        
                if subj == 3 && (trial == 2 || trial == 4) % special case where force data is not trustworthy, do not store, ok to do calc's for debugging
                    TrialData(n).Results.Force = nan;
                    
                    % x/ML dir
                    TrialData(n).Results.cx_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqx_torso = nan; %
                    TrialData(n).Results.px_torso = nan; %
                    TrialData(n).Results.cx_torso_st = nan; % regression coeff's
                    TrialData(n).Results.rsqx_torso_st = nan; %
                    TrialData(n).Results.cx_lag_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqx_lag_torso = nan; %
                    TrialData(n).Results.px_lag_torso = nan; %
                    
                    TrialData(n).Results.cx_IP = nan; % regression coeff's
                    TrialData(n).Results.rsqx_IP = nan; %
                    TrialData(n).Results.px_IP = nan; %
                    
                    TrialData(n).Results.IntCumWork = nan;
                    
                    TrialData(n).Results.corrPowerFresLagX = nan;
                    TrialData(n).Results.corrPFxLagpval = nan;
                    TrialData(n).Results.corrPFxpval = nan;
                    TrialData(n).Results.FresLagX = nan;
                    
                    % z/Vert dir
                    TrialData(n).Results.cz_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqz_torso = nan; %
                    TrialData(n).Results.pz_torso = nan; %
                    TrialData(n).Results.cz_torso_st = nan; % regression coeff's
                    TrialData(n).Results.rsqz_torso_st = nan; %
                    
                    % 2D vector xz dir
                    TrialData(n).Results.FxzNorm = nan;
                    TrialData(n).Results.theta = nan;
                    TrialData(n).Results.cxz_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqxz_torso = nan; %
                    TrialData(n).Results.pxz_torso = nan; %
                    TrialData(n).Results.cxz_torso_st = nan; % regression coeff's
                    TrialData(n).Results.rsqxz_torso_st = nan; %
                    
                    % Power
                    TrialData(n).Results.IntPower = nan; % x dir
                    TrialData(n).Results.IntPowerXZ = nan; % 2D
                    TrialData(n).Results.POBPower = nan; % x dir
                    TrialData(n).Results.POBPowerXZ = nan; % 2D
                    
                    % Angular - some of these metrics are calculated/stored
                    % earlier in code. Now overwrite bad trials' results
                    % with nan's
                    TrialData(n).Results.TyTorso = nan;
                    TrialData(n).Results.TyGd = nan;
                    TrialData(n).Results.Tg = nan;
                    TrialData(n).Results.TyExt = nan;
                    TrialData(n).Results.TyInt = nan;
                    TrialData(n).Results.TyIntXhZh = nan;
                    TrialData(n).Results.TyIntXlZh = nan;
                    TrialData(n).Results.TyIntXhZl = nan;
                    TrialData(n).Results.TyIntXlZl = nan;
                    TrialData(n).Results.TyFxh = nan;
                    TrialData(n).Results.TyFxl = nan;
                    TrialData(n).Results.TyFzh = nan;
                    TrialData(n).Results.TyFzl = nan;
                    TrialData(n).Results.c_ang_torso = nan; % regression coeff's
                    TrialData(n).Results.rsq_ang_torso = nan; % rsquare
                    TrialData(n).Results.p_ang_torso = nan; % p-value
                    TrialData(n).Results.c_angFx_torso = nan; % regression coeff's
                    TrialData(n).Results.rsq_angFx_torso = nan; % rsquare
                    TrialData(n).Results.p_angFx_torso = nan; % p-value
                    TrialData(n).Results.c_angFz_torso = nan; % regression coeff's
                    TrialData(n).Results.rsq_angFz_torso = nan; % rsquare
                    TrialData(n).Results.p_angFz_torso = nan; % p-value
                    TrialData(n).Results.c_angInt_torso = nan; % regression coeff's
                    TrialData(n).Results.rsq_angInt_torso = nan; % rsquare
                    TrialData(n).Results.p_angInt_torso = nan; % p-value
                    TrialData(n).Results.angP = nan;
                    TrialData(n).Results.lagTyAngTorso = nan;
                    TrialData(n).Results.xcorrTyAngTorso = nan;
                    TrialData(n).Results.lagTyWTorso = nan;
                    TrialData(n).Results.xcorrTyWTorso = nan;
                    TrialData(n).Results.lagTyAlphaTorso = nan;
                    TrialData(n).Results.xcorrTyAlphaTorso = nan;
                    
                else % good force data
                
                    if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                        % Arm vec's
                        TrialData(n).Results.AssistArm = armAsst;
                        TrialData(n).Results.AssistArmVel = vArmAsst;
                        TrialData(n).Results.AssistArmAcc = aArmAsst;
                        % Store POB arm length data
                        TrialData(n).Results.POBArm = armPOB;
                        TrialData(n).Results.POBArmVel = vArmPOB;
                        TrialData(n).Results.POBArmAcc = aArmPOB;
                        % Save interaction point state
                        TrialData(n).Results.IntPt = lFin;
%                         TrialData(n).Results.IntPt = rFin;
                        TrialData(n).Results.IntPtVel = vFin;
                        TrialData(n).Results.IntPtAcc = aFin;
                        
                        % Save info on fit force to pos, vel, acc of torso marker
                        % x/ML dir
                        TrialData(n).Results.cx_torso = cx_torso; % regression coeff's
                        TrialData(n).Results.rsqx_torso = rsqx_torso; %
                        TrialData(n).Results.px_torso = px_torso; %
                        TrialData(n).Results.cx_torso_st = cx_torso_st; % regression coeff's
                        TrialData(n).Results.rsqx_torso_st = rsqx_torso_st; %
%                         TrialData(n).Results.cx_lag_torso = cx_lag_torso; % regression coeff's
%                         TrialData(n).Results.rsqx_lag_torso = rsqx_lag_torso; %
%                         TrialData(n).Results.px_lag_torso = px_lag_torso; %
                        % Save info on fit force to IP
%                         TrialData(n).Results.cx_IP = cx_IP; % regression coeff's
%                         TrialData(n).Results.rsqx_IP = rsqx_IP; %
%                         TrialData(n).Results.px_IP = px_IP; %
%                         TrialData(n).Results.cy_torso = cy_torso; % regression coeff's
%                         TrialData(n).Results.rsqy_torso = rsqy_torso; % rsquare
%                         TrialData(n).Results.py_torso = py_torso; % p-value
                        % z/Vert dir
                        TrialData(n).Results.cz_torso = cz_torso; % regression coeff's
                        TrialData(n).Results.rsqz_torso = rsqz_torso; % rsquare
                        TrialData(n).Results.pz_torso = pz_torso; % p-value
                        
                        % 2D vector xz dir
                        TrialData(n).Results.FxzNorm = Fmag;
                        TrialData(n).Results.theta = theta;
                        TrialData(n).Results.cxz_torso = cxz_torso; % regression coeff's
                        TrialData(n).Results.rsqxz_torso = rsqxz_torso; %
                        TrialData(n).Results.pxz_torso = pxz_torso; %
                        TrialData(n).Results.cxz_torso_st = cxz_torso_st; % regression coeff's
                        TrialData(n).Results.rsqxz_torso_st = rsqxz_torso_st; %
                        
                        % Angular
                        TrialData(n).Results.c_ang_torso = c_ang_torso; % regression coeff's
                        TrialData(n).Results.rsq_ang_torso = rsq_ang_torso; % rsquare
                        TrialData(n).Results.p_ang_torso = p_ang_torso; % p-value
                        TrialData(n).Results.c_angFx_torso = c_angFx_torso; % regression coeff's
                        TrialData(n).Results.rsq_angFx_torso = rsq_angFx_torso; % rsquare
                        TrialData(n).Results.p_angFx_torso = p_angFx_torso; % p-value
                        TrialData(n).Results.c_angFz_torso = c_angFz_torso; % regression coeff's
                        TrialData(n).Results.rsq_angFz_torso = rsq_angFz_torso; % rsquare
                        TrialData(n).Results.p_angFz_torso = p_angFz_torso; % p-value
                        TrialData(n).Results.c_angInt_torso = c_angInt_torso; % regression coeff's
                        TrialData(n).Results.rsq_angInt_torso = rsq_angInt_torso; % rsquare
                        TrialData(n).Results.p_angInt_torso = p_angInt_torso; % p-value
                                                
                        % Power
                        TrialData(n).Results.IntPower = IP_power; % x dir
                        TrialData(n).Results.IntPowerXZ = IP_power_XZ; % 2D
                        TrialData(n).Results.POBPower = POB_power; % x dir
                        TrialData(n).Results.POBPowerXZ = POB_power_XZ; % 2D
    %                     TrialData(n).Results.IPpowerX_lo = IP_powerX_lo;
    %                     TrialData(n).Results.IPpowerX_hi = IP_powerX_hi;
                                        % Angular power
                                        
                        TrialData(n).Results.IntCumWork = IP_cumWork;
                    end
                end
                
                
                clear IP_power IP_cumWork 
            end
            
            %% Check plots
            % Plot to check that start and end times are correct, also check height of midline/reference
            if plotCheck == 1 
                t = 0:(length(torso(:,2))-1); 
                t = t./sample_rate;
                if any(strcmpi(TrialData(n).Info.Condition,{'Solo Ground','Assist Ground'}))
                    plotind = plotind + 1;
                    if plotind == 1
                        f1 = figure;
                    end
                    set(0, 'currentfigure', f1); 
                    subplot(numrows,numcols,plotind)
                    yyaxis left
                    plot(t,torso(:,1)-TrialData(n).Results.midline,'g'), hold on;
                    ylim([-0.2 0.2]);hline(0,'k--');
%                         plot(2:length(t),vTorso(:,2),'g'), hold on;
                    if ismember(plotind,[1 6 11 16])
%                         ylabel('AP disp (m)');
                        ylabel('Torso ML pos (m)')
                    end
                    yyaxis right
                elseif any(strcmpi(TrialData(n).Info.Condition,{'Solo Beam','Assist Beam','Assist beam'}))
                    plot2ind = plot2ind + 1;
                    if plot2ind == 1
                        f2 = figure;
                    end
                    set(0, 'currentfigure', f2);
                    subplot(numrows,numcols,plot2ind)
                end

                plot(t,LHEE(:,3)/1000,'b-',t,RHEE(:,3)/1000,'r-'),hold on;
                ylim([0 0.4]);
%                     plot(1:length(t),LHEE(:,3)/1000,'b-',1:length(t),RHEE(:,3)/1000,'r-'),hold on;
                vline([t(start_idx) t(stop_idx)],'k-');  
%                     vline([start_idx  stop_idx],'k-');  
                hline(TrialData(n).Results.ht,'k-');
                if ismember(plotind,numcols*(1:4))
                    ylabel('Vert pos (m)');
                end
                if plotind == 1 && gcf == f1
                    legend('Torso','Lheel','Rheel','orientation','horizontal')
                elseif plot2ind == 1 && gcf == f2
                    legend('Lheel','Rheel','orientation','horizontal')
                end
                xlabel('Time (s)');
                if (plotind == 1 && gcf == f1) || (plot2ind == 1 && gcf == f2)
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);   
            elseif plotCheck == 2 && any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam','Assist Beam','Assist beam'})) % plot distance calculated from heel vert algo vs. manually entered distance traveled. 
                % Must plot beam distance here in order to compare Luke's method vs.
                % algo method. Plot solo in diff color from assist so see if learning
                % trend across all trials or just solo trials
                plotind = plotind + 1; % not for subplots, just for labels
                hold on;
                tn = str2num(TrialData(n).Info.Trial(end-1:end));
                if strcmpi(TrialData(n).Info.Condition,'Solo Beam')
                    plot(tn,distA(n),'bx'),...
                    plot(tn,distL(n)/1000,'bo')
                elseif strcmpi(TrialData(n).Info.Condition,'Assist Beam','Assist beam')
                    plot(tn,distA(n),'kx'),...
                    plot(tn,distL(n)/1000,'ko')
                end
                if plotind == 1
                    legend('algo','entered'); ylabel('Distance completed (m)'),xlabel('Trial');
                    set(gcf,'outerposition',[672   662   576   389]);
                    titlename = sprintf('HHI%i',subj); title(titlename);
                end
            elseif plotCheck == 3 
                t = 0:(length(torso(:,2))-1); 
                t = t./sample_rate;
                if any(strcmpi(TrialData(n).Info.Condition,{'Solo Ground','Assist Ground'}))
                    plotind = plotind + 1;
                    if plotind == 1
                        f1 = figure;
                    end
                    set(0, 'currentfigure', f1); 
                    subplot(numrows,numcols,plotind)
                    plot(t,COM(:,1)-TrialData(n).Results.midline,'b-'),hold on; % Plot COM for overgd to see if different from torso and which better to use for midline
                elseif any(strcmpi(TrialData(n).Info.Condition,{'Solo Beam','Assist Beam','Assist beam'}))
                    plot2ind = plot2ind + 1;
                    if plot2ind == 1
                        f2 = figure;
                    end
                    set(0, 'currentfigure', f2);
                    subplot(numrows,numcols,plot2ind)
                end
                plot(t,torso(:,1)-TrialData(n).Results.midline,'k-')
                % Must plot lines after full plot so scaled correctly
                ylim([-0.5 0.5]);
                hline(0,'k--');
                vline([t(start_idx) t(stop_idx)],'k-');  

                if ismember(plotind,numcols*(1:4))
                    ylabel('ML pos (m)');
                end
                xlabel('Time (s)');
                if (plotind == 1 && gcf == f1) || (plot2ind == 1 && gcf == f2)
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                    if (plotind == 1 && gcf == f1)
                        legend('COM','Torso');
                    end
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);  
            elseif plotCheck == 4 && any(strcmpi(TrialData(n).Info.Condition,{'Assist Beam','Assist beam','Assist Ground'})) % Plot force and torque on L axis and ma ground (beam or floor) on R only for assisted trials
                t = 0:(length(Force(:,2))-1); 
                t = t./sample_rate;
                clear Fxy
                for i = 1:length(Force)
                    Fxy(i) = norm(Force(i,[1 3]));
                end

                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                yyaxis left
                plot(t,Fxy,t,Fxy'.*TrialData(n).Results.maGd,t,TrialData(n).Results.TyGd), hold on;
                if subj == 5 && trial == 2
                    ylim([-300 300]);
                else
                    ylim([-20 20]);
                end
                if ismember(plotind,[1 6 11 16 21])
                    ylabel('|F| or Torque (N or Nm)')
                end
                hline(0,'k-');

                yyaxis right
                plot(t,TrialData(n).Results.maGd)
                ylim([-1 1]);
                if ismember(plotind,[1 6 11 16 21]+4)
                    ylabel('Moment arm gd (m)')
                end

%                     vline([t(start_idx) t(stop_idx)],'k-');  
%                     vline([start_idx  stop_idx],'k-');  

                if plotind == 1 
                    legend('|Fxy|','|Fxy|*m.a.','Ty','m.a.','orientation','horizontal')
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);

                xlabel('Time (s)');
            elseif plotCheck == 5 && any(strcmpi(TrialData(n).Info.Condition,{'Solo Beam','Solo beam','Assist Beam','Assist beam'})) % Plot w (torso angular speed) to check angular momentum
                t = 0:(length(TrialData(n).Results.wTorso)-1); 
                t = t./sample_rate;

                plotind = plotind + 1;
                subplot(numrows,numcols,plotind)
                if any(strcmpi(TrialData(n).Info.Condition,{'Assist Beam','Assist beam'}))
                    yyaxis left
                end
                plot(t,TrialData(n).Results.wTorso),hold on;    
                ylim([-3 3]*.1);
                if ismember(plotind,[1 6 11 16 21])
                    ylabel('Torso ang vel (rad/s)')
                end

                if any(strcmpi(TrialData(n).Info.Condition,{'Assist Beam','Assist beam'}))
                    yyaxis right
                    plot(t,TrialData(n).Results.TyGd)
                    ylabel('Torque (Nm)'),
                    if subj == 5 && trial == 2
                        ylim([-200 200]);
                    else
                        ylim([-20 20]);
                    end
                end

%                     vline([t(start_idx) t(stop_idx)],'k-');  
%                     vline([start_idx  stop_idx],'k-');  
                hline(0,'k-');

                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);

                xlabel('Time (s)');
            elseif plotCheck == 6 & any(strcmpi(TrialData(n).Info.Condition,{'Assist Ground','Assist Beam','Assist beam'})) & ~isnan(TrialData(n).Results.Force)
                firstPlot = firstPlot + 1;
                if firstPlot == 1
                    fx = figure; fz = figure; ty_fx = figure; ty_fz = figure; ty_int = figure;
                end
                plotind = plotind + 1;
                % Fx
                figure(fx);
                subplot(numrows,numcols,plotind)
                histogram(TrialData(n).Results.Force(:,1),'Normalization','probability');
                set(gca,'tickdir','out','box','off');
                if mod(plotind,numcols)==1
                    ylabel('Probability')
                end
                if plotind > 15
                    xlabel('Fx (N)');
                end
                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);
                % Fz
                figure(fz);
                subplot(numrows,numcols,plotind)
                histogram(TrialData(n).Results.Force(:,3),'Normalization','probability');
                set(gca,'tickdir','out','box','off');
                if mod(plotind,numcols)==1
                    ylabel('Probability')
                end
                if plotind > 15
                    xlabel('Fz (N)');
                end
                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);
                % Ty_Fx
                figure(ty_fx);
                subplot(numrows,numcols,plotind)
                histogram(TrialData(n).Results.TyFx,'Normalization','probability');
                set(gca,'tickdir','out','box','off');
                if mod(plotind,numcols)==1
                    ylabel('Probability')
                end
                if plotind > 15
                    xlabel('TyFx (Nm)');
                end
                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);
                % Ty_Fz
                figure(ty_fz);
                subplot(numrows,numcols,plotind)
                histogram(TrialData(n).Results.TyFz,'Normalization','probability');
                set(gca,'tickdir','out','box','off');
                if mod(plotind,numcols)==1
                    ylabel('Probability')
                end
                if plotind > 15
                    xlabel('TyFz (Nm)');
                end
                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);
                % Ty_Int
                figure(ty_int);
                subplot(numrows,numcols,plotind)
                histogram(TrialData(n).Results.TyInt,'Normalization','probability');
                set(gca,'tickdir','out','box','off');
                if mod(plotind,numcols)==1
                    ylabel('Probability')
                end
                if plotind > 15
                    xlabel('TyInt (Nm)');
                end
                if plotind == 1 
                    titlename = sprintf('HHI%i %s %s',subj,TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                else
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                end
                title(titlename);
            end
                
            % RESULTS FOR ALL TRIALS
            TrialData(n).Results.time = time; % only the analysis window
            TrialData(n).Results.startIdx = start_idx;
            TrialData(n).Results.stopIdx = stop_idx;
        end
    end
end

function [start,stop] = checkArm(Markers,start,stop)

% Check and make sure both 'arm' markers are present for the entire
% start/stop window
RSHO = Markers.RSHO(start:stop,1);
RFIN = Markers.RFIN(start:stop,1);

% Finger or Shoulder must be at least 100mm (10cm) away from the origin of
% the Vicon Coordinate System
idxSHO = abs(RSHO) > 100;
idxFIN = abs(RFIN) > 100;

startSHO = find(idxSHO,1,'first');
stopSHO = find(idxSHO,1,'last');

startFIN = find(idxFIN,1,'first');
stopFIN = find(idxFIN,1,'last');
% Recalculate start and stop times
start = max([start,startSHO+start-1,startFIN+start-1]);
stop = min([stop,stopSHO+start-1,stopFIN+start-1]);
% Adjust for ringing from filtering.
if startSHO ~= 1 || startFIN ~= 1
    start = start+5;
end
if stopSHO ~= length(RSHO) || stopFIN ~= length(RSHO)
   stop = stop - 5; 
end
end