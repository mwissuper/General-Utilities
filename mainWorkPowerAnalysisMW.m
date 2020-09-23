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
    numrows = 4; numcols = 5;
    plotind = 0;
    
    if subj == 4 
        TrialData(16) = []; % bad processing in Nexus
    elseif subj == 7 
        TrialData(51) = []; % trial to repeat trial 23 bc no gait mat in trial 23. But doesn't matter for analysis, trial 23 looks fine.
    elseif subj == 14
        TrialData(26) = []; % notes say to skip this and replace with trial 51. trial looks fine in Nexus though...
    end
    % Number of trials in the dataset:
    num_trials = length(TrialData);
    distA = nan*ones(size(num_trials)); distL = nan*ones(size(num_trials));
    
    %% Find midline of beam for reference. Do it here so can use for
    % later reference. Cycle through all beam-walking trials to
    % get array to plot. Calculate mean as midline of beam.
    temp.beamMidline = [];
    plotind = 0;
    for n = 1:num_trials
        trial = str2num(TrialData(n).Info.Trial(end-1:end));
        if ~isempty(strfind(TrialData(n).Info.Condition,'Beam')) || ~isempty(strfind(TrialData(n).Info.Condition,'beam')) % if beam-walking trial
            if ~isempty(strfind(TrialData(n).Info.Condition,'Assist')) % if 2 person marker data set
                Markers = medfiltFields(TrialData(n).Markers.POB,1);
            else
                Markers = medfiltFields(TrialData(n).Markers,1);
            end
            % Find start of trial as when foot on ground at beg of trial
            % gets to its first peak
            LHEE = Markers.LHEE;
            RHEE = Markers.RHEE;
            vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
            %% Get POB torso state values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
            if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                torsoY = Markers.C7(:,2)./1000;
            else
                torsoY = Markers.CLAV(:,2)./1000;
            end
            [start_idx,stop_idx] = getHHIAnalysisWindow_MW(Markers,vLHEEZfilt,1,torsoY);
            % Special exceptions based on plots and find point by eye.
            if subj == 3 
                if trial == 27 || trial == 29
                    start_idx = 500;
                elseif trial == 49
                    start_idx = 430;
                end
            elseif subj == 4
                if trial == 1
                    start_idx = 500;
                elseif trial == 7
                    start_idx = 600;
                elseif trial == 22
                    start_idx = 712;
                elseif trial == 26
                    start_idx = 620;
                elseif trial == 29 || trial == 35
                    start_idx = 500;
                elseif trial == 46
                    start_idx = 600;
                end
            elseif subj == 5
                if trial == 3
                    start_idx = 930;
                elseif trial == 7
                    start_idx = 812;
                elseif trial == 18
                    start_idx = 365;
                elseif trial == 20
                    start_idx = 497;
                elseif trial == 23
                    start_idx = 200;
                elseif trial == 29
                    start_idx = 314;
                elseif trial == 37
                    start_idx = 306;
                elseif trial == 49
                    start_idx = 270;
                end
            elseif subj == 8
                if trial == 28
                    start_idx = 240;
                end
            elseif subj == 9
                if trial == 5
                    start_idx = 290;
                elseif trial == 6
                    start_idx = 600;
                elseif trial == 21
                    start_idx = 261;
                elseif trial == 43
                    start_idx = 240;
                end
            elseif subj == 11
                if trial == 12
                    start_idx = 240;
                end
            elseif subj == 12
                if trial == 9 || trial == 24 || trial == 25 || trial == 35
                    start_idx = 500;
                elseif trial == 43
                    start_idx = 470;
                end
            elseif subj == 13
                if trial == 1
                    start_idx = 840;
                elseif trial == 16
                    start_idx = 500;
                elseif trial == 17
                    start_idx = 430;
                elseif trial == 23
                    start_idx = 350;
                end
            end
            
%             % Plot to check RHEE marker at start_idx looks reasonable for beam
%             % midline (i.e. no issues with gaps)
%             plotind = plotind + 1;
%             subplot(4,5,plotind)
%             plot(RHEE(:,1)),hold on,plot(start_idx,RHEE(start_idx,1),'x');
%             title(TrialData(n).Info.Trial); ylabel('RHEE ML pos (m)');
            
            temp.beamMidline = [temp.beamMidline RHEE(start_idx,1)];
            
            TrialData(n).Results.beamMidline = temp.beamMidline;
        end 
    end
    beamMidline = nanmean(temp.beamMidline)/1000;
    beamMidlineSD = nanstd(temp.beamMidline)/1000;
    
    clear Markers
    temp = [];

    %% Trial Analysis Main Loop
    for n = 1:num_trials
        disp(['Processing ',TrialData(n).Info.Trial]);
        trial = str2num(TrialData(n).Info.Trial(end-1:end));
        % Check for trials to process
        if any(strcmpi(TrialData(n).Info.Condition, to_process))
            %% -- Compute start and stop time of the trial. --- %%
            clear start_idx stop_idx
            % Median filter all of the marker data to remove jumps. Then
            % calculate the start and stop indices for HHI analysis
            max_distance = TrialData(n).Info.Distance_Traveled.*25.4; %convert from in to mm
            if strcmpi(TrialData(n).Info.Condition, 'Assist Solo') % Assist Solo
                Markers = TrialData(n).Markers;
                if isfield(Markers,'AP') == 1 % Somet trials denoted assistant by AP while others did not
                    Markers.AP = medfiltFields(Markers.AP,1);
                    LHEE = Markers.AP.LHEE;
                    RHEE = Markers.AP.RHEE;
                    vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
                    %% Get AP torso values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
                    if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                        torsoY = Markers.AP.C7(:,2)./1000;
                    else
                        torsoY = Markers.AP.CLAV(:,2)./1000;
                    end
                    vTorsoY = filtfilt(sos, g, diff(torsoY)); % Filter before use for algo to find analysis window
                    % Luke's method using max_distance
                    [start_L,stop_L] = getHHIAnalysisWindow(Markers.AP,max_distance);
                    % Check the assister's arm
                    [start_L,stop_L] = checkArm(Markers.AP,start_L,stop_L);

                    % New method
                    [start_idx,stop_idx] = getHHIAnalysisWindow_MW(Markers.AP,vLHEEZfilt,0,vTorsoY);
                else
                    Markers = medfiltFields(TrialData(n).Markers,1);
                    LHEE = Markers.LHEE;
                    RHEE = Markers.RHEE;
                    vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
                    %% Get AP torso values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
                    if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                        torsoY = Markers.C7(:,2)./1000;
                    else
                        torsoY = Markers.CLAV(:,2)./1000;
                    end
                    vTorsoY = filtfilt(sos, g, diff(torsoY)); % Filter before use for algo to find analysis window
                    % Luke's method
                    [start_L,stop_L] = getHHIAnalysisWindow(Markers,max_distance);
                    % Check the assister's arm
                    [start_L,stop_L] = checkArm(Markers,start_L,stop_L);

                    [start_idx,stop_idx] = getHHIAnalysisWindow_MW(Markers,vLHEEZfilt,0,vTorsoY);
                end
            elseif any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist Ground','Assist beam'}))
                % Assist Ground, Assist Beam
                Markers = TrialData(n).Markers;
                Markers.AP = medfiltFields(Markers.AP,1);
                Markers.POB = medfiltFields(Markers.POB,1);
                LHEE = Markers.POB.LHEE;
                RHEE = Markers.POB.RHEE;
                vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
                %% Get POB torso values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
                if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                    torsoY = Markers.POB.C7(:,2)./1000;
                else
                    torsoY = Markers.POB.CLAV(:,2)./1000;
                end
                vTorsoY = filtfilt(sos, g, diff(torsoY)); % Filter before use for algo to find analysis window
                % Luke's method
                [start_L,stop_L] = getHHIAnalysisWindow(Markers.POB,max_distance);
                % Check the assister's arm
%                 [start_L,stop_L] = checkArm(Markers.AP,start_L,stop_L);
                
                % My method
                if strcmpi(TrialData(n).Info.Condition, 'Assist Ground')
                    [start_idx,stop_idx] = getHHIAnalysisWindow_MW(Markers.POB,vLHEEZfilt,0,vTorsoY); % use max distance method and not heel vertical height relative to beam when no beam-walking
                else
                    vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
                    [start_idx,stop_idx,beamHt] = getHHIAnalysisWindow_MW(Markers.POB,vLHEEZfilt,1);
                end
            else % Solo
                % Solo Beam, Solo Ground
                Markers = medfiltFields(TrialData(n).Markers,1);
                LHEE = Markers.LHEE;
                RHEE = Markers.RHEE;
                vLHEEZfilt = filtfilt(sos, g, diff(LHEE(:,3))); 
                %% Get POB torso values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
                if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                    torsoY = Markers.C7(:,2)./1000;
                else
                    torsoY = Markers.CLAV(:,2)./1000;
                end
                % Luke's method
                [start_L,stop_L] = getHHIAnalysisWindow(Markers,max_distance);
                if strcmpi(TrialData(n).Info.Condition, 'Assist Ground')
%                     start_idx = start_L; stop_idx = stop_L; % use max distance method and not heel vertical height relative to beam when no beam-walking
                    [start_idx,stop_idx] = getHHIAnalysisWindow_MW(Markers,vLHEEZfilt,0,vTorsoY);
                else
                    [start_idx,stop_idx,beamHt] = getHHIAnalysisWindow_MW(Markers,vLHEEZfilt,1,torsoY);
                end
            end
            
            %% Special exception trials that can't be fixed algorithmically
            % found index by eye in Nexus
            if subj == 3 
                if trial == 2
                    start_idx = start_idx + 20; % transient at beginning due to incorrect gapfill
                    stop_idx = 1357;
                elseif trial == 3
                    stop_idx = 419; % took 2 tries in this trial, just count first one
                elseif trial == 4
                    stop_idx = 1447;
                elseif trial == 13 
                    stop_idx = 441;
                elseif trial == 20
                    stop_idx = 829;
                elseif trial == 29
                    stop_idx = 1329;
                elseif trial == 36
                    stop_idx = 1143;
                elseif trial == 37
                    stop_idx = 442;
                elseif trial == 42
                    stop_idx = 1376;
                elseif trial == 46
                    stop_idx = 1204;
                end
            elseif subj == 4
                if trial == 1
                    stop_idx = 1350;
                elseif trial == 3
                    stop_idx = 3792;
                    start_idx = start_idx + 50; % transient at beginning due to incorrect gapfill
                elseif trial == 7
                    stop_idx = 2972;
                elseif trial == 9
                    stop_idx = 1773;
                elseif trial == 11
                    stop_idx = 3908;
                elseif trial == 17 
                    stop_idx = 1960;
                elseif trial == 22
                    stop_idx = 2224;
                elseif trial == 24
                    stop_idx = 2014;
                elseif trial == 26
                    stop_idx = 1776;
                elseif trial == 35
                    stop_idx = 1321;
                elseif trial == 39
                    stop_idx = 1413;
                elseif trial == 40
                    stop_idx = 1388;
                elseif trial == 43
                    stop_idx = 1520;
                elseif trial == 45
                    stop_idx = 1105;
                elseif trial == 46
                    stop_idx = 1261;
                end
            elseif subj == 5
                if trial == 3
                    stop_idx = 779;
                elseif trial == 4
                    stop_idx = 1537;
                elseif trial == 7
                    stop_idx = 885;
                elseif trial == 13
                    stop_idx = 685;
                elseif trial == 17
                    stop_idx = 1179;
                elseif trial == 18
                    stop_idx = 960;
                elseif trial == 21
                    stop_idx = 467;
                elseif trial == 23
                    stop_idx = 985;
                elseif trial == 29
                    stop_idx = 824;
                elseif trial == 35
                    stop_idx = 815;
                elseif trial == 36
                    stop_idx = 796;
                elseif trial == 37
                    stop_idx = 943;
                elseif trial == 42
                    stop_idx = 534;
                elseif trial == 44
                    stop_idx = 489;
                elseif trial == 46
                    stop_idx = 842;
                elseif trial == 49
                    stop_idx = 661;
                end
            elseif subj == 7
                if trial == 3
                    stop_idx = 1436;
                elseif trial == 7
                    stop_idx = 1176;
                elseif trial == 20
                    stop_idx = 1139;
                elseif trial == 21
                    stop_idx = 997;
                elseif trial == 27
                    stop_idx = 1146;
                elseif trial == 37
                    stop_idx = 963;
                end
            elseif subj == 8
                if trial == 5
                    stop_idx = 1042;
                elseif trial == 8
                    stop_idx = 987;
                elseif trial == 9
                    stop_idx = 807;
                elseif trial == 10
                    stop_idx = 712;
                elseif trial == 11
                    stop_idx = 953;
                elseif trial == 13
                    stop_idx = 807;
                elseif trial == 15
                    stop_idx = 810;
                elseif trial == 21
                    stop_idx = 774;
                elseif trial == 23
                    stop_idx = 838;
                elseif trial == 24
                    stop_idx = 764;
                elseif trial == 28
                    stop_idx = 772;
                elseif trial == 36
                    stop_idx = 822;
                elseif trial == 39
                    stop_idx = 592;
                elseif trial == 42
                    stop_idx = 634;
                elseif trial == 46
                    stop_idx = 757;
                elseif trial == 47
                    stop_idx = 814;
                elseif trial == 48
                    stop_idx = 693;
                end
            elseif subj == 9
                if trial == 1
                    stop_idx = 941;
                elseif trial == 5
                    stop_idx = 1659;
                elseif trial == 6
                    stop_idx = 1564;
                elseif trial == 17
                    stop_idx = 1265;
                elseif trial == 20
                    stop_idx = 1313;
                elseif trial == 22
                    stop_idx = 1396;
                elseif trial == 28
                    stop_idx = 1258;
                elseif trial == 33
                    stop_idx = 1173;
                elseif trial == 39
                    stop_idx = 1202;
                elseif trial == 46
                    stop_idx = 1212;
                elseif trial == 49
                    stop_idx = 879;
                elseif trial == 50
                    stop_idx = 1093;
                end
            elseif subj == 10
                if trial == 3
                    stop_idx = 1554;
                elseif trial == 11
                    stop_idx = 807;
                elseif trial == 23
                    stop_idx = 1116;
                elseif trial == 24
                    stop_idx = 948;
                elseif trial == 45
                    stop_idx = 1028;
                elseif trial == 47
                    stop_idx = 885;
                end
            elseif subj == 13
                if trial == 9
                    stop_idx = 549;
                elseif trial == 10
                    stop_idx = 1322;
                elseif trial == 34
                    stop_idx = 558;
                end
            elseif subj == 14
                if trial == 8
                    stop_idx = 437;
                end
            end

            %% Construct the time axis
            sample_rate = TrialData(n).Markers.samplerate;
            time = 0:(stop_idx - start_idx);
            time = time./sample_rate;
            % Save Lheel pos used for finding window
            LHEEZ = LHEE(:,3);
            RHEEZ = RHEE(:,3); 
            % Calculations depending on force data, 2 person marker set
            if ~isempty(strfind(TrialData(n).Info.Condition, 'Assist'))
                %% Window and Transpose Force and Marker Data %% 
                % Now that we have the start and stop indices, we can get the
                % force data and the marker data for the time when the subject
                % is walking.
                
                % The Markers we need are the Right Shoulder annd Finger of
                % the Assistance Provider (m)
                if isfield(Markers,'AP') == 1                   
                    RSHO = Markers.AP.RSHO(start_idx:stop_idx, :)./1000;
                    RFIN = Markers.AP.RFIN(start_idx:stop_idx, :)./1000;
                else
                    RSHO = Markers.RSHO(start_idx:stop_idx, :)./1000;
                    RFIN = Markers.RFIN(start_idx:stop_idx, :)./1000;
                end
                
                if isfield(Markers,'POB') == 1
                    LSHO = Markers.POB.LSHO(start_idx:stop_idx, :)./1000;
                    if isfield(Markers.POB,'LFIN') == 1
                        LFIN = Markers.POB.LFIN(start_idx:stop_idx, :)./1000; 
                    else
                        LFIN = nan;
                    end
                end

                % Get the sampling rate for the forces                
                force_rate = TrialData(n).Forces.samplerate;
                
                clear Force Torque
                % Collect the forces and torques into the initialized arrays
                for k = 1:3 % for each direction
                    %Downsample the forces and torques to be at the same
                    %sampling rate as the marker data
                    F = resample(TrialData(n).Forces.(force_names{k}), sample_rate, force_rate);
                    T = resample(TrialData(n).Forces.(torque_names{k}), sample_rate, force_rate);
                    % Median filter forces to reduce noise
                    F = medfilt1(F);
                    T = medfilt1(T);
                    % Limit the forces and torques to the same analysis window
                    % as the marker data
                    Force(:,k) = F(start_idx:stop_idx);
                    Torque(:,k) = T(start_idx:stop_idx);
                end
                
                % Checked for all trials with avail video assist beam
                % Fx volt is negative when in shear away from center of
                % sensor, so must negate sign before processing. Also
                % negate Fy volt bc that is consistent with setup for
                % BL experiment fall 2019. Do this before recoverForces
                % projects force signals into lab CS. 3/30/20 MW.
                Force(:,1) = -Force(:,1);
                Force(:,2) = -Force(:,2);
                
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
                
                % Need to subtract out force of gravity on FH for vertical
                % dir and subtract out force of accelerating mass of FH in
                % each direction (?)
                
%% --------------- WORK and POWER ANALYSIS ------------------------------%%
                % Now, filter the force vector to obtain smooth
                % results. 
                clear rFin rFilt vFin vFilt RASIFilty armAsst armPOB POBpower
                for k = 1:3
                    Force(:,k) = filtfilt(sos, g, Force(:,k));
                    rFilt(:,k) = filtfilt(sos, g, RFIN(:,k)); 
                end
                rFin = RFIN; % Don't double-filter!
 
                % Look at each direction separately. Must filter before
                % take deriv. Again, here the signs are reversed for force
                % bc it's force on POB
                vFin = diff(rFilt).*sample_rate;
                for i = 1:3
                    intPt_power(:,i) = Force(2:end,i).*vFin(:,i); % P < 0 for motor because F > 0 for tension and v > 0 for move to the right
                    if i == 1 % x axis, check worst case drift effect of 1.5N
                        IP_powerX_hi(1,:) = (Force(2:end,i)+1.5).*vFin(:,i);
                        IP_powerX_lo(1,:) = (Force(2:end,i)-1.5).*vFin(:,i);
                    end
                    intPt_cumWork(:,i) = cumsum(intPt_power(:,i))./sample_rate;
                    intPt_netWork_norm(i) = sum(intPt_power(:,i))/sample_rate/(time(end)-time(1)); % One number characterizing whole trial overall, normalize by time length of trial
                end
                
                %% Calc xcorr Fint to variables. Fit Fint to variables
                
                for i = 1:3
                    vFilt(:,i) = filtfilt(sos, g, vFin(:,i));
                end
                aFin = diff(vFilt).*sample_rate;               
                
                if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                    clear torso torsoFilt vTorso vTorsoFilt aTorso
                    %% Get POB torso state values. Use clav but C7 if needed. Abs coordinate system, not subtract out midline.
                    if (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                        torso = Markers.POB.C7(start_idx:stop_idx,:)./1000;
                    else
                        torso = Markers.POB.CLAV(start_idx:stop_idx,:)./1000;
                    end
                    for k = 1:3
                        torsoFilt(:,k) = filtfilt(sos, g, torso(:,k)); % (m) % Use only for taking deriv.
                    end

                    % Store torso pos for later analysis/plotting
                    TrialData(n).Results.torso = torso; % (m)
                    vTorso = diff(torsoFilt).*sample_rate; 
                    % Filter vel before take derivative
                    for k = 1:3
                        vTorsoFilt(:,k) = filtfilt(sos, g, vTorso(:,k));
                    end
                    aTorso = diff(vTorsoFilt).*sample_rate;
                    % Store torso derivatives for later analysis/plotting
                    TrialData(n).Results.vTorso = vTorso;
                    TrialData(n).Results.aTorso = aTorso;

                    %% Calculate xcorr for a variety of signals to compare to LT work and to adjust regression if needed
                    % for xcorr(x,y), y lags x if lag < 0;
                    
                    % FIP to torso disp
                    % Must remove means before doing xcorr. Not matter if
                    % remove midline of beam or walking path correctly for
                    % xcorr, removing mean no matter what
                    temp.Fx = Force(:,1) - nanmean(Force(:,1));
                    temp.swayX = torso(:,1) - nanmean(torso(:,1)); % This is not displacement! This is position with mean removed!
                    [r, lags] = xcorr(temp.Fx,temp.swayX,sample_rate,'normalized'); % Constrain window for xcorr to +/- 1s
                    ind = find(abs(r) == max(abs(r))); % Look mag of xcorr
                    if length(ind) > 1 % Not sure what to do here, just flag it for now
                        msg = sprintf('More than one max for xcorr!');
                    end
                    TrialData(n).Results.lagFIPTorsoX = lags(ind)/sample_rate;
                    TrialData(n).Results.xcorrFIPTorsoX = r(ind);
                    clear r lags ind
                  
                    % FIP to POB vTorso in ML dir only                    
                    [r, lags] = xcorr(temp.Fx(2:end),vTorso(:,1),sample_rate,'normalized'); % Constrain window for xcorr to +/- 1s
                    ind = find(abs(r) == max(abs(r)));
                    TrialData(n).Results.lagFIPvTorsoX = lags(ind)/sample_rate;
                    TrialData(n).Results.xcorrFIPvTorsoX = r(ind);
                    clear r lags ind
                    
                    % vIP to POB vTorso in ML dir only. Use vIn that's not
                    % double-filtered.                 
                    [r, lags] = xcorr(vFin(:,1),vTorso(:,1),sample_rate,'normalized'); % Constrain window for xcorr to +/- 1s
                    ind = find(abs(r) == max(abs(r)));
                    TrialData(n).Results.lagvIPvTorsoX = lags(ind)/sample_rate;
                    TrialData(n).Results.xcorrvIPvTorsoX = r(ind);
                    clear r lags ind
                    temp = [];
                    
                    %% Fit Fint to torso state var's 

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
                    % away from IP. For fit model, -Fx/y/z is consistent with a
                    % positive spring or damper restoring person to midline
                    % or starting position
                    clear m;
                    lag = floor(.157*sample_rate); % calculated from xcorr analysis 
                    i = 1; % Just do for x dir. Should write code to use beamHt reference if want to fit model in vertical/z dir. Not make sense to fit model to AP/y dir.
                    ref = beamMidline; 
                    % No lag
                    m = [aTorso(:,i) vTorso(2:end,i) torso(3:end,i)-ref ones(size(aTorso(:,i)))]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                    % Lag of 157ms
                    mlag = [aTorso(lag:end,i) vTorso(lag+1:end,i) torso(lag+2:end,i)-ref ones(size(aTorso(lag:end,i)))]; 

                    [cx_torso,rsqx_torso,px_torso] = regressIter(Force(3:end,i),m); % Fx > 0 if tension. If tension and displacement > 0 (POB move to R), should get positive spring
                    [cx_lag_torso,rsqx_lag_torso,px_lag_torso] = regressIter(Force(lag+2:end,i),mlag); 
                    % Test correlations
                    % no lag
                    TrialData(n).Results.FresX = Force(3:end,i) - m*cx_torso;  % Calculate corr power to residual. However, oppositte force sign conventions for power vs. force for model, so just know that -rho values mean correlated.
                    [temp, TrialData(n).Results.corrPFxpval] = corr(TrialData(n).Results.FresX,intPt_power(2:end,i));
                    if TrialData(n).Results.corrPFxpval < 0.05
                        TrialData(n).Results.corrPowerFresX = temp;
                    else
                        TrialData(n).Results.corrPowerFresX = nan;
                    end
                    % lag
                    TrialData(n).Results.FresLagX = Force(lag+2:end,i) - mlag*cx_lag_torso;  % Calculate corr power to residual. However, oppositte force sign conventions for power vs. force for model, so just know that -rho values mean correlated.
                    [temp, TrialData(n).Results.corrPFxLagpval] = corr(TrialData(n).Results.FresLagX,intPt_power(lag+1:end,i));
                    if TrialData(n).Results.corrPFxLagpval < 0.05
                        TrialData(n).Results.corrPowerFresLagX = temp;
                    else
                        TrialData(n).Results.corrPowerFresLagX = nan;
                    end               
  
                    clear temp;
                    
                    %% Fit Fint to IP state var's 

                    % Regress per direction. Regress
                    % iteratively, throwing away any variables with coeff's
                    % whose CI's include zero. For ML disp, look at IP
                    % relative to a midline. Assume model acts to
                    % restore IP to mean IP pos in ML dir. Must reverse
                    % sign of force if want positive damper and spring coeff
                    % to mean passive elements. 
                    % Up to this point, sign convention for all force
                    % directions is positive when tension/shear for POB
                    % away from IP. For fit model, -Fx/y/z is consistent with a
                    % positive spring or damper restoring person to midline
                    % or starting position
                    clear m;
                    i = 1; % Just do for x dir. Should write code to use IP height reference if want to fit model in vertical/z dir. Not make sense to fit model to AP/y dir.
                    ref = nanmean(rFin(:,1)); 
                    % No lag
                    m = [aFin(:,i) vFin(2:end,i) rFin(3:end,i)-ref ones(size(aFin(:,i)))]; % Want displacement, not abs pos for fitting. 

                    [cx_IP,rsqx_IP,px_IP] = regressIter(Force(3:end,i),m); % Fx > 0 if tension. If tension and displacement > 0 (POB move to R), should get positive spring          
  
                    clear temp;
                    
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
                
%% --------------- STORE THE RESULTS ----------------------------------- %%
                % For assisted trials, we store the forces, the arm vector,
                % and the work, power, and alignment
                        
                if subj == 3 && (trial == 2 || trial == 4) % special case where force data is not trustworthy, do not store, ok to do calc's for debugging
                    TrialData(n).Results.Forces = nan;
                    
                    TrialData(n).Results.cx_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqx_torso = nan; %
                    TrialData(n).Results.px_torso = nan; %
                    TrialData(n).Results.cx_lag_torso = nan; % regression coeff's
                    TrialData(n).Results.rsqx_lag_torso = nan; %
                    TrialData(n).Results.px_lag_torso = nan; %
                    
                    TrialData(n).Results.cx_IP = nan; % regression coeff's
                    TrialData(n).Results.rsqx_IP = nan; %
                    TrialData(n).Results.px_IP = nan; %
                    
                    TrialData(n).Results.IntPower = nan;
                    TrialData(n).Results.IntCumWork = nan;
                    
                    TrialData(n).Results.corrPowerFresLagX = nan;
                    TrialData(n).Results.corrPFxLagpval = nan;
                    TrialData(n).Results.corrPFxpval = nan;
                    TrialData(n).Results.FresLagX = nan;
                else
                    TrialData(n).Results.Forces = Force;
                
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
                        TrialData(n).Results.IntPt = rFin;
                        TrialData(n).Results.IntPtVel = vFin;
                        TrialData(n).Results.IntPtAcc = aFin;
                        % Save info on fit force to pos, vel, acc of torso marker
                        TrialData(n).Results.cx_torso = cx_torso; % regression coeff's
                        TrialData(n).Results.rsqx_torso = rsqx_torso; %
                        TrialData(n).Results.px_torso = px_torso; %
                        TrialData(n).Results.cx_lag_torso = cx_lag_torso; % regression coeff's
                        TrialData(n).Results.rsqx_lag_torso = rsqx_lag_torso; %
                        TrialData(n).Results.px_lag_torso = px_lag_torso; %
%                         TrialData(n).Results.cy_torso = cy_torso; % regression coeff's
%                         TrialData(n).Results.rsqy_torso = rsqy_torso; % rsquare
%                         TrialData(n).Results.py_torso = py_torso; % p-value
%                         TrialData(n).Results.cz_torso = cz_torso; % regression coeff's
%                         TrialData(n).Results.rsqz_torso = rsqz_torso; % rsquare
%                         TrialData(n).Results.pz_torso = pz_torso; % p-value
                        % Save info on fit force to IP
                        TrialData(n).Results.cx_IP = cx_IP; % regression coeff's
                        TrialData(n).Results.rsqx_IP = rsqx_IP; %
                        TrialData(n).Results.px_IP = px_IP; %
                    end

                    TrialData(n).Results.IntPower = intPt_power;
                    TrialData(n).Results.IPpowerX_lo = IP_powerX_lo;
                    TrialData(n).Results.IPpowerX_hi = IP_powerX_hi;
                    TrialData(n).Results.IntCumWork = intPt_cumWork;
                end
                
                
                clear intPt_power intPt_cumWork 
            else % one person data set, no force data
                clear torso torsoFilt vTorso vTorsoFilt aTorso
                %% Get POB torso state values. Use clav or C7 if needed
                if (subj == 4 && trial == 7) || (subj == 5 && (trial == 21 || trial == 31 || trial == 42)) || (subj == 7 && (trial == 13 || trial == 44)) || (subj == 14 && trial == 8) 
                    if isfield(Markers,'AP') % Assist Solo cond's
                        torso = Markers.AP.C7(start_idx:stop_idx,:)./1000;
                    else
                        torso = Markers.C7(start_idx:stop_idx,:)./1000;
                    end
                else
                    if isfield(Markers,'AP') % Assist Solo cond's
                        torso = Markers.AP.CLAV(start_idx:stop_idx,:)./1000;
                    else
                        torso = Markers.CLAV(start_idx:stop_idx,:)./1000;
                    end
                end
                
                for k = 1:3
                    torsoFilt(:,k) = filtfilt(sos, g, torso(:,k)); % (m) % Use only for taking deriv.
                end

                % Store torso pos for later analysis/plotting
                TrialData(n).Results.torso = torso; % (m)
                vTorso = diff(torsoFilt).*sample_rate; 
                % Filter vel before take derivative
                for k = 1:3
                    vTorsoFilt(:,k) = filtfilt(sos, g, vTorso(:,k));
                end
                aTorso = diff(vTorsoFilt).*sample_rate;
                % Store torso derivatives for later analysis/plotting
                TrialData(n).Results.vTorso = vTorso;
                TrialData(n).Results.aTorso = aTorso;
            end
            %% Calculate kinematic metrics like distance completed and sway
            % Use the medfiltered data from Markers struct
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam','Assist beam','Solo Beam', 'Assist Ground', 'Solo Ground'}))
                sway = torso(:,1);
                if any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam', 'Solo Ground'}))                   
                    temp.LPSI = Markers.LPSI(start_idx:stop_idx,1)./1000;
                    temp.RPSI = Markers.RPSI(start_idx:stop_idx,1)./1000;
                    COMLatPos = mean([temp.LPSI'; temp.RPSI'],1);
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
                    temp.LPSI = Markers.POB.LPSI(start_idx:stop_idx,1)./1000;
                    temp.RPSI = Markers.POB.RPSI(start_idx:stop_idx,1)./1000;
                    COMLatPos = mean([temp.LPSI; temp.RPSI]);
                    % Linear correlation between lateral COM sway and lateral
                    % forces and correlation bewteen lateral torso sway and
                    % vertical forces. Need to remove nan's before do corr.
                    ind = ~isnan(COMLatPos);
                    [rho, TrialData(n).Results.pFlatCOM] = corr(Force(ind,1),COMLatPos(ind));
                    if TrialData(n).Results.pFlatCOM < 0.05
                        TrialData(n).Results.rFlatCOM = rho;
                    else
                        TrialData(n).Results.rFlatCOM = nan;
                    end
                    ind = ~isnan(sway);
                    [rho, TrialData(n).Results.pFsway] = corr(Force(ind,1),sway(ind));
                    if TrialData(n).Results.pFsway < 0.05
                        TrialData(n).Results.rFsway = rho;
                    else
                        TrialData(n).Results.rFsway = nan;
                    end

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
                
                if ~isempty(strfind(TrialData(n).Info.Condition,'beam'))
                    midline = beamMidline;
                else
                    midline = nanmean(sway);
                end
                sway = sway - midline; % same as torso pos with midline removed
                TrialData(n).Results.beamerSway = sway;
                TrialData(n).Results.beamerCOMSway = COMLatPos;
                
                % Distance is calculated as the difference between the
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
                %% Plot to check that start and end times are correct
                % Plot whole trial but xlim to start and stop times with a
                % little before and after. Plot hline for beamht
                if plotCheck == 1 && any(strcmpi(TrialData(n).Info.Condition, {'Solo Ground','Assist Ground'}))%  {'Solo Beam','Assist Beam','Assist beam'} plot to check distance completed on beam looks reasonable. Plot AP and ML pos of torso relative to beam center
                    plotind = plotind + 1;
                    subplot(numrows,numcols,plotind)
                    t = 0:(length(torso(:,2))-1); % time starts and ends with start_idx and stop_idx
                    t = t./sample_rate;
                    yyaxis left
%                     plot(t,(left_heel(:,2)-left_heel(start_idx,2))/1000,t,(right_heel(:,2)-left_heel(start_idx,2))/1000),hold on; % subtract init value so easy to compare to plots of dispA and dispL
                    plot(t,(torso(:,2)-torso(start_idx,2))/1000),hold on; % subtract init value so easy to compare to plots of dispA and dispL
                    if any(plotind,[1 6 11 16])
                        ylabel('AP disp (m)');
                    end
%                     xlim([t(start_idx)-1 t(stop_idx)+1]);
%                     vline([t(start_idx) t(stop_idx)],'k-');
                    yyaxis right
                    plot(t,LHEEZ(start_idx:stop_idx)/1000,t,RHEEZ(start_idx:stop_idx)/1000);
%                     hline(beamHt/1000,'k-'); % convert to m to plot. Only
%                     plot for beam-walking trials
                    if any(plotind,numcols*(1:4))
                        ylabel('Vert pos (m)');
                    end
                    if plotind == 1
                        legend('Torso','Lheel','Rheel','orientation','horizontal')
                    end
                    xlabel('Time (s)');
                    if plotind == 1
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
                elseif plotCheck == 3 && any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam','Assist Beam','Assist beam'}))  % Check that beam midline and SD look reasonable for each trial's sway
                    plotind = plotind + 1;
                    subplot(numrows,numcols,plotind)
                    t = 0:(length(torso(:,2))-1); % time starts and ends with start_idx and stop_idx
                    t = t./sample_rate;
                    plot(t,torso(:,1)),hold on;
                    hline(beamMidline,'k');
                    hline(beamMidline + [-beamMidlineSD beamMidlineSD],'k--');
                    xlabel('Time (s)'),ylabel('Torso pos (m)');
                    if plotind == 1
                        titlename = sprintf('HHI%i %s',subj,TrialData(n).Info.Trial); title(titlename);
                    else
                        titlename = sprintf('HHI%s',TrialData(n).Info.Trial); title(titlename);
                    end
                elseif plotCheck == 4 && any(strcmpi(TrialData(n).Info.Condition, {'Solo Ground','Assist Ground'})) % Check time window of analysis by plot sway and vertical displacement (visually vert disp is a strong tell of when they stopped walking)
                    plotind = plotind + 1;
                    subplot(numrows,numcols,plotind)
                    t = 0:(length(torso(:,2))-1); % time starts and ends with start_idx and stop_idx
                    t = t./sample_rate;
                    yyaxis left
                    plot(t,torso(:,1)),hold on;
                    xlabel('Time (s)'),ylabel('Torso ML pos (m)');
                    yyaxis right
                    plot(t(2:end),vTorsoFilt(:,2)),ylabel('Torso AP vel (m)'),ylim([0 0.005]);
                    if plotind == 1
                        titlename = sprintf('HHI%i %s',subj,TrialData(n).Info.Trial); title(titlename);
                    else
                        titlename = sprintf('HHI%s',TrialData(n).Info.Trial); title(titlename);
                    end
                end
                
                %% We also store the average speed of the trial, and a
                % logical value indicating if the trial was completed.
                TrialData(n).Results.avgSpeed = TrialData(n).Results.totalDistance./time(end); % time has been trimmed to cover start to stop only
                TrialData(n).Results.completed = (TrialData(n).Info.Distance_Traveled == 144);
            elseif strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                TrialData(n).Results.totalDistance = nan;
                distA(n) = nan; distL(n) = nan;
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