function [TrialData] = mainWorkPowerAnalysisMW(TrialData,subj)
    %MAINWORKPOWERANALYSIS: Processing for HHI Experimental Data collected after
    %August 2017
    %
    %   TrialData = mainWorkPowerAnalysis(TrialData) performs a series of analyses on
    %   the data stored in TrialData. mainWorkPowerAnalysis loops through the data
    %   stored in TrialData and processes trials of type 'Assist
    %   Ground', 'Assist Beam', or 'Assist Solo.' The analysis includes:
    %       1) Recovering the interaction forces in Vicon Coordinates
    %       2) Calculating Work and Power
    %   mainWorkPowerAnalysis also processes trials of type 'Solo Beam'; however,
    %   these trials do not require a work analysis
    %
    %   mainWorkPowerAnalysis adds a field called Results to TrialData. For Assist
    %   Beam and Assist Solo Trials, the fields of Results are:
    %       AssistArm: The arm length vector for the assistance provider
    %       (mm)
    %       AssistWork: the work performed by the Assister (J)
    %       AssistPower: the power transferred by the Assister (W)
    %       time: the time vector for the forces and arm vector (s)
    %       PowerAlignment: the alignment between the resultant force and
    %       the assister arm velocity vector. Alignment is in [-1,1] and is
    %       alignment = cos(theta) = dot(F,v)./norm(F)*norm(v) (nu)
    %       startidx: the collection index for which analysis starts (nu)
    %       stopidx: the collection index for which the analysis ends (nu)
    %   For Assist Beam and Solo Beam Trials, the Results field contains:
    %       beamerSway: the CLAV lateral displacement referenced to the start
    %       CLAV position at the start of the trial (mm)
    %       totalDistance: the total distance the beamer traveled along the
    %       beam(mm)
    %       avgSpeed: the average speed of the person on the beam
    %       throughout the trial (mm/s)
    %       completed: a boolean indicating if the person on the beam
    %       completed the walk (equivalent to traveling 12 ft along the
    %       beam)
    %   Note: Force calibration and conversion is performed in a separate
    %   file recoverForces.m. Adjustments to recoverForces.m may be needed
    %   as the force sensor is recalibrated. Additionally, recoverForces.m
    %   does not take into account the weight of the force handle, which
    %   may be of imporance later on.
    
    %   Luke Drnach
    %   October 30, 2017
    
    %   Refactored on December 4, 2018. Individual routines were separated
    %   into functions for code clarity and potential future re-use.
    
    %   Sometimes people turned their heads near end of trial and seemed to
    %   be getting ready to change direction, don't want that to affect
    %   metrics so reduce max_distance traveled by fixed amount and then 
    %   checked individual trial plots of power peaks to see not drop to ~0
    %   at end of trial. MW 11/1/19
    distOffset = 250; % mm
    
    %   MW made lots of changes to remove filtering on marker data and do
    %   minimal filtering on force data. The marker data has likely been
    %   filtered a great deal already through procbatch.
    
    %% Initialization
    % Number of trials in the dataset:
    num_trials = length(TrialData);
    % Trial Types to Process
    to_process = {'Assist Beam','Assist Ground','Solo Beam','Solo Ground','Assist Solo'};
    % Assist Trial Names
    assisted = {'Assist Beam','Assist Ground','Assist Solo'};
    % Field names for force and torque
    force_names = {'Fx','Fy','Fz'};
    torque_names = {'Mx','My','Mz'};
    
    % Mass of FT sensor
    mFT = 4/9.81; % (kg)
    
    %% Trial Analysis Main Loop
    for n = 1:num_trials
        disp(['Processing ',TrialData(n).Info.Trial]);
        % Check for trials to process
        if any(strcmpi(TrialData(n).Info.Condition, to_process))
            % Get sampling freq, set up filters
            % Butterworth 3rd order lowpass filter, cutoff 20Hz for markers 60hz
            % for force
            fs = TrialData(n).Markers.samplerate;
            [Bm,Am] = butter(3,20/(fs/2)); % Use only to filter vel before diff to get acc
            fsa = TrialData(n).Forces.samplerate;
            [Bf,Af] = butter(3,60/(fsa/2)); % Using this before downsample 
            
            %% -- Compute start and stop time of the trial. --- %%
            
            % Median filter all of the marker data to remove jumps. Then
            % calculate the start and stop indices for HHI analysis
            max_distance = TrialData(n).Info.Distance_Traveled*25.4 - distOffset; % Convert in to mm, subtract out extra in case of head turn
            if strcmpi(TrialData(n).Info.Condition, 'Assist Solo') % Assist Solo
                Markers = TrialData(n).Markers;
                if isfield(Markers,'AP') == 1
                    Markers.AP = medfiltFields(Markers.AP,1);
                    [start_idx,stop_idx] = getHHIAnalysisWindow(Markers.AP,max_distance);
                else
                    Markers = medfiltFields(TrialData(n).Markers,1);
                    [start_idx,stop_idx] = getHHIAnalysisWindow(Markers,max_distance);
                end
            elseif any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Assist Ground'}))
                % Assist Ground, Assist Beam
                Markers = TrialData(n).Markers;
                Markers.AP = medfiltFields(Markers.AP,1);
                Markers.POB = medfiltFields(Markers.POB,1);
                [start_idx,stop_idx] = getHHIAnalysisWindow(Markers.POB,max_distance);
                % Check the assister's arm
                [start_idx,stop_idx] = checkArm(Markers.AP,start_idx,stop_idx);
            else
                % Solo Beam, Solo Ground
                Markers = medfiltFields(TrialData(n).Markers,1);
                [start_idx,stop_idx] = getHHIAnalysisWindow(Markers,max_distance);
            end
            % Special exception participants where we need to set stop_idx
            % (end of trial) manually since code fails
            if subj == 4
                if strcmp(TrialData(n).Info.Trial,'Trial16')
                    stop_idx = start_idx + 750;
                end
            elseif subj == 10
                if strcmp(TrialData(n).Info.Trial,'Trial46')
                    stop_idx = start_idx + 450;
                end
            end

            %Construct the time axis
            time = 0:(stop_idx - start_idx);
            time = time./fs;
            if any(strcmpi(TrialData(n).Info.Condition, assisted))
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
               
                clear Force Torque F T
                % Collect the forces and torques into the initialized arrays
                for k = 1:3 % for each direction
                    %Filter forces and torques before downsampling
                    F = filtfilt(Bf,Af,TrialData(n).Forces.(force_names{k}));
                    T = filtfilt(Bf,Af,TrialData(n).Forces.(torque_names{k}));
                    %Downsample the forces and torques to be at the same
                    %sampling rate as the marker data
                    F = resample(F, fs, fsa);
                    T = resample(T, fs, fsa);
                    % Limit the forces and torques to the same analysis window
                    % as the marker data
                    Force(:,k) = F(start_idx:stop_idx);
                    Torque(:,k) = T(start_idx:stop_idx);
                end
                
                % Checked for all trials with avail video assist beam
                % Fx volt is negative when in shear away from center of
                % sensor, so must negate sign before processing. Also
                % negative Fy volt bc that is consistent with setup for
                % BL experiment fall 2019. Do this before recoverForces
                % projects force signals into lab CS. 3/30/20 MW.
                Force(:,1) = -Force(:,1);
                Force(:,2) = -Force(:,2);
                
                if isfield(Markers,'FH')
                    % Median filter force handle MARKERS to remove jumps.
                    % These markers were not previously medfiltered
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
                % Filter the force vector 
                clear rFin rFilt vFin vFilt aFin clav clavFilt vClav vClavFilt aClav RASIFilty armAsst armPOB POBpower

                % Look at each direction separately. The signs are reversed for force
                % bc it's force on POB
                vFin = diff(RFIN).*fs;
                for i = 1:3
                    intPt_power(:,i) = -Force(2:end,i).*vFin(:,i); % If Fx > 0, IP in shear with partners' force vectors away from each other. if vx > 0, IP is towards the POB on the right. By this convention, Px > 0 means brake. We want Px > 0 to mean motor
                    intPt_cumWork(:,i) = cumsum(intPt_power(:,i))./fs;
                    intPt_netWork_norm(i) = sum(intPt_power(:,i))/fs/(time(end)-time(1)); % One number characterizing whole trial overall, normalize by time length of trial
                end
                
                %% Fit Fint to int pt disp, vel, acc data for each direction - old analysis. Not really useful. Fit Fint to clavicle is more interesting.
                % Filter velocity before diff to get acc for intPt 
                for i = 1:3
                    vFilt(:,i) = filtfilt(Bm, Am, vFin(:,i));
                end
                aFin = diff(vFilt).*fs;               
                
                if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                    %% Get POB clavicle state values
                    clav = Markers.POB.CLAV(start_idx:stop_idx,:)./1000;

                    % Store clavicle pos for later analysis/plotting
                    TrialData(n).Results.CLAV = clav; % (m)
                    vClav = diff(clav).*fs; 
                    % Filter vel before take derivative
                    for k = 1:3
                        vClavFilt(:,k) = filtfilt(Bm, Am, vClav(:,k));
                    end
                    aClav = diff(vClavFilt).*fs;
                    % Store clavicle derivatives for later analysis/plotting
                    TrialData(n).Results.vCLAV = vClav;
                    TrialData(n).Results.aCLAV = aClav;

%                     %% Look at POB as rigid body and take power as force at
%                     % hands multiplied by clav velocity
%                     for i = 1:3
%                         POBpower(:,i) = Force(2:end,i).*vClav(:,i);
%                     end
%                     TrialData(n).Results.POBpower = POBpower; % Store results
%                     
                    %% Fit Fint to clavicle state var's 

                    % Regress per direction. Regress
                    % iteratively, throwing away any variables with coeff's
                    % whose CI's include zero. For ML disp, look at clav
                    % relative to beam midline, assumed we can get it from
                    % starting ML pos of LHEE (code to get analysis window
                    % always looks at LHEE for start). Assume model acts to
                    % restore clav to middle of beam in ML dir and to clav
                    % height and sagittal pos at beg of trial. Must reverse
                    % sign of force if we're looking at force felt by POB.
                    % Then we get correct signs for damper and spring (i.e.
                    % Up to this point, sign convention for all force
                    % directions is positive when tension or shear for POB
                    % away from IP. For fit model, -Fx/y/z is consistent with a
                    % positive spring or damper restoring person to midline
                    % or starting position
                    % of beam. The force we're looking at here is not the
                    % same as the force up to this point [which was force
                    % on IP]. Now we're looking at force on POB.)
                    clear m;
                    for i = 1:3
                        [trough, itrough] = findpeaks(-Markers.POB.LHEE(start_idx:end,3),'MinPeakProminence',range(Markers.POB.LHEE(start_idx:end,3)/5)); % find troughs after first peak
                        ind = itrough(1) + start_idx - 1;
                        if i == 1 % Find beam midline by pos LHEE at first local min after start_idx (start_idx is first vert peak LHEE) 
                            ref = Markers.POB.LHEE(ind,1)/1000; 
                            TrialData(n).Results.beamMidline = ref; % Store this value for plotting later. Just x component
                        else
                            ref = clav(ind,i); % pos of clav when first step on to beam
                        end
                        m = [aClav(:,i) vClav(2:end,i) clav(3:end,i)-ref ones(size(aClav(:,i)))]; % Want displacement, not abs pos for fitting. For now, look at changes in pos of interaction point
                        if i == 1
                            [cx_clav,rsqx_clav,px_clav] = regressIter(Force(3:end,i),m); % Fx > 0 if tension. If tension and displacement > 0 (POB move to R), should get positive spring
                            TrialData(n).Results.FresX = Force(3:end,i) - m*cx_clav;  % Calculate corr power to residual. However, oppositte force sign conventions for power vs. force for model, so just know that -rho values mean correlated.
                            [temp, TrialData(n).Results.corrPFxpval] = corr(TrialData(n).Results.FresX,intPt_power(2:end,i));
                            if TrialData(n).Results.corrPFxpval < 0.05
                                TrialData(n).Results.corrPowerFresX = temp;
                            else
                                TrialData(n).Results.corrPowerFresX = nan;
                            end
                        elseif i == 2
                            [cy_clav,rsqy_clav,py_clav] = regressIter(Force(3:end,i),m); 
                        else
                            [cz_clav,rsqz_clav,pz_clav] = regressIter(Force(3:end,i),m); 
                        end               
                    end   
                    clear temp;

                    %% Calculate xcorr Fint to clavicle disp just to compare to LT work
                    
                    % Must remove means before doing xcorr
                    temp.Fx = Force(:,1) - nanmean(Force(:,1));
    %                 temp.swayX = clavFilt(1,:) - nanmean(clavFilt(1,:));
                    temp.swayX = clav(:,1) - nanmean(clav(:,1)); % This is not displacement! This is position with mean removed!
                    [r, lags] = xcorr(temp.Fx,temp.swayX,fs,'normalized'); % Constrain window for xcorr to +/- 1s
                    ind = find(abs(r) == max(abs(r)));
                    if length(ind) > 1 % Not sure what to do here, just flag it for now
                        ind = ind
                    end
                    lagX = lags(ind)/fs;
                    xcorrX = r(ind);

                    temp.Fz = Force(:,3) - nanmean(Force(:,3));
                    temp.swayZ = clav(:,3) - nanmean(clav(:,3));
                    [r, lags] = xcorr(temp.Fx,temp.swayZ,fs,'normalized');
                    ind = find(abs(r) == max(abs(r)));
                    if length(ind) > 1 % Not sure what to do here, just flag it for now
                        ind = ind
                    end
                    lagZ = lags(ind)/fs;
                    xcorrZ = max(r);

                    temp = [];

                    %% Get POB RASI marker to calculate AP velocity of whole body 
                    % Need to filter first before taking derivative!
                    RASIFilty = filtfilt(Bm, Am, Markers.POB.RASI(start_idx:stop_idx,2))./1000; % (m)
                    % Store pos for later analysis/plotting
                    TrialData(n).Results.RASI = Markers.POB.RASI(start_idx:stop_idx,2)./1000;

                    % Store velocity for later analysis/plotting
                    TrialData(n).Results.vyRASI = diff(RASIFilty).*fs;
                    
                    %% Get Assistant's arm len state 
                end
                
%% --------------- STORE THE RESULTS ----------------------------------- %%
                % For assisted trials, we store the forces, the arm vector,
                % and the work, power, and alignment
                TrialData(n).Results.Forces = Force;
                
                if ~strcmpi(TrialData(n).Info.Condition, 'Assist Solo')
                    
                    % For assisted trials, we store the forces, the arm vector,
                    % and the work, power, and alignment

                    % Store POB arm length data

                    
%                     % Save info on fit force to pos, vel, acc of rFin marker
%                     TrialData(n).Results.cx = cx; % regression coeff's
%                     TrialData(n).Results.rsqx = rsqx; %
%                     TrialData(n).Results.px = px; %
%     %                 TrialData(n).Results.cintx = cintx; % confidence intervals on regression coeff's
%                     TrialData(n).Results.cy = cy; % regression coeff's
%                     TrialData(n).Results.rsqy = rsqy; % rsquare
%                     TrialData(n).Results.py = py; % p-value
%     %                 TrialData(n).Results.cinty = cinty; % confidence intervals on regression coeff's
%                     TrialData(n).Results.cz = cz; % regression coeff's
%                     TrialData(n).Results.rsqz = rsqz; % rsquare
%                     TrialData(n).Results.pz = pz; % p-value
%     %                 TrialData(n).Results.cintz = cintz; % confidence intervals on regression coeff's

                    % Save info on fit force to pos, vel, acc of clavicle marker
                    TrialData(n).Results.cx_clav = cx_clav; % regression coeff's
                    TrialData(n).Results.rsqx_clav = rsqx_clav; %
                    TrialData(n).Results.px_clav = px_clav; %
                    TrialData(n).Results.cy_clav = cy_clav; % regression coeff's
                    TrialData(n).Results.rsqy_clav = rsqy_clav; % rsquare
                    TrialData(n).Results.py_clav = py_clav; % p-value
                    TrialData(n).Results.cz_clav = cz_clav; % regression coeff's
                    TrialData(n).Results.rsqz_clav = rsqz_clav; % rsquare
                    TrialData(n).Results.pz_clav = pz_clav; % p-value

                    % Save xcorr data
                    TrialData(n).Results.xcorrX = xcorrX; 
                    TrialData(n).Results.lagX = lagX; 
                    TrialData(n).Results.xcorrZ = xcorrZ;
                    TrialData(n).Results.lagZ = lagZ; 
                end
                
                % Save interaction point state and work/power
                TrialData(n).Results.IntPt = RFIN;
                TrialData(n).Results.IntPtVel = vFin;
                TrialData(n).Results.IntPtAcc = aFin;
%                 TrialData(n).Results.IntCumWork_tot = intPt_cumWork_tot;
%                 TrialData(n).Results.IntPower_tot = intPt_power_tot;
%                 TrialData(n).Results.IntAlignment = intPt_alignment;
                TrialData(n).Results.IntPower = intPt_power;
                TrialData(n).Results.IntCumWork = intPt_cumWork;
                
                clear intPt_power intPt_cumWork xcorrX lagX xcorrZ lagZ
            end
            %% Sway and segment angle 
            % For beam walking trials, we also store the POB mediolateral
            % sway and the distance traveled along the beam as well as
            % force and sway correlations
            % Also collect ML sway as calculated by midpoing PSIS markers
            % Use the medfiltered data from Markers struct
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Beam', 'Solo Beam', 'Assist Ground', 'Solo Ground'}))
                clear clavFilt vClavFilt aClav
                if any(strcmpi(TrialData(n).Info.Condition, {'Solo Beam', 'Solo Ground'}))
                    % Store clavicle pos for later analysis/plotting
                    sway = Markers.CLAV(start_idx:stop_idx,1)./1000; % (m)
                    clavY = Markers.CLAV(:,2)./1000; % (m)
                    temp.LPSI = Markers.LPSI(start_idx:stop_idx,1)./1000;
                    temp.RPSI = Markers.RPSI(start_idx:stop_idx,1)./1000;
                    COMsway = mean([temp.LPSI'; temp.RPSI'],1);
                    TrialData(n).Results.corrLat = nan;
                    TrialData(n).Results.corrVert = nan;
                    % Also collect thorax (approx with
                    % projection of line between two thorax markers) ang
                    % rotation for POB. Calculate leg segment obliq by
                    % 1) project all malleolus markers and PSIS markers
                    % into frontal plane, 2) take midpoint 2 malleoli, 3)
                    % take midpoint PSIS markers, 4) take vector between
                    % the two points and get angle in frontal plane
                    TrialData(n).Results.pelvicObliq = TrialData(n).Markers.LPelvisAngles(start_idx:stop_idx,1)./1000; % L and R pelvis angles are exact mirror images of each other. Pelvic data not really necessary, just look at thorax and leg segments for now     
                    temp.u = Markers.C7(start_idx:stop_idx,:) - Markers.CLAV(start_idx:stop_idx,:); % Vector from clavicle to C7 as approx thorax vertical
                    temp.u = temp.u./1000;
                    % Leg segment angle (only care x and z components since
                    % will look at frontal plane angle
                    temp.mANK(:,1) = mean([Markers.LANK(start_idx:stop_idx,1) Markers.RANK(start_idx:stop_idx,1)],2)./1000;
                    temp.mANK(:,3) = mean([Markers.LANK(start_idx:stop_idx,3) Markers.RANK(start_idx:stop_idx,3)],2)./1000;
                    temp.mPSIS(:,1) = mean([Markers.LPSI(start_idx:stop_idx,1) Markers.RPSI(start_idx:stop_idx,1)],2)./1000;
                    temp.mPSIS(:,3) = mean([Markers.LPSI(start_idx:stop_idx,3) Markers.RPSI(start_idx:stop_idx,3)],2)./1000;
                    temp.legVec = temp.mPSIS-temp.mANK;
                    v = [0 0 1]; % Vertical line
%                     for i = 1:length(temp.u)
%                         TrialData(n).Results.thoraxObliq(i) = getAngVec(temp.u(i,:),v,2);
%                         TrialData(n).Results.legObliq(i) = getAngVec(temp.legVec(i,:),v,2);
%                     end  
                    temp = [];
                else % Force data exists, two people marker set
                    sway = Markers.POB.CLAV(start_idx:stop_idx,1)./1000;
                    clavY = Markers.POB.CLAV(:,2)./1000;
                    temp.LPSI = Markers.POB.LPSI(start_idx:stop_idx,1)./1000;
                    temp.RPSI = Markers.POB.RPSI(start_idx:stop_idx,1)./1000;
                    COMsway = mean([temp.LPSI; temp.RPSI]);
                    % Linear correlation between lateral COM sway and lateral
                    % forces and correlation bewteen lateral Clav sway and
                    % vertical forces. Need to remove nan's before do corr.
                    ind = ~isnan(COMsway);
                    [rho, TrialData(n).Results.pFlatCOM] = corr(Force(ind,1),COMsway(ind));
                    if TrialData(n).Results.pFlatCOM < 0.05
                        TrialData(n).Results.rFlatCOM = rho;
                    else
                        TrialData(n).Results.rFlatCOM = nan;
                    end
                    ind = ~isnan(sway);
                    [rho, TrialData(n).Results.pFvertClav] = corr(Force(ind,3),sway(ind));
                    if TrialData(n).Results.pFvertClav < 0.05
                        TrialData(n).Results.rFvertClav = rho;
                    else
                        TrialData(n).Results.rFvertClav = nan;
                    end
                    % Look at other two combo's force and sways just to
                    % check corr
                    ind = ~isnan(COMsway);
                    [rho, TrialData(n).Results.pFvertCOM] = corr(Force(ind,3),COMsway(ind));
                    if TrialData(n).Results.pFvertCOM < 0.05
                        TrialData(n).Results.rFvertCOM = rho;
                    else
                        TrialData(n).Results.rFvertCOM = nan;
                    end
                    ind = ~isnan(sway);
                    [rho, TrialData(n).Results.pFlatClav] = corr(Force(ind,1),sway(ind));
                    if TrialData(n).Results.pFlatClav < 0.05
                        TrialData(n).Results.rFlatClav = rho;
                    else
                        TrialData(n).Results.rFlatClav = nan;
                    end
                    % Also collect thorax (approx with
                    % projection of line between two thorax markers) ang
                    % rotation for POB
                    if isfield(Markers.POB,'LPelvisAngles')
                        TrialData(n).Results.pelvicObliq = TrialData(n).Markers.POB.LPelvisAngles(start_idx:stop_idx,1); % L and R pelvis angles are exact mirror images of each other  
                    else
                        TrialData(n).Results.pelvicObliq = nan;
                    end
                    temp.u = Markers.POB.C7(start_idx:stop_idx,:) - Markers.POB.CLAV(start_idx:stop_idx,:); % Vector from clavicle to C7 as approx thorax vertical
                    temp.u = temp.u./1000;
                    % Leg segment angle (only care x and z components since
                    % will look at frontal plane angle
                    temp.mANK(:,1) = mean([Markers.POB.LANK(start_idx:stop_idx,1) Markers.POB.RANK(start_idx:stop_idx,1)],2)./1000;
                    temp.mANK(:,3) = mean([Markers.POB.LANK(start_idx:stop_idx,3) Markers.POB.RANK(start_idx:stop_idx,3)],2)./1000;
                    temp.mPSIS(:,1) = mean([Markers.POB.LPSI(start_idx:stop_idx,1) Markers.POB.RPSI(start_idx:stop_idx,1)],2)./1000;
                    temp.mPSIS(:,3) = mean([Markers.POB.LPSI(start_idx:stop_idx,3) Markers.POB.RPSI(start_idx:stop_idx,3)],2)./1000;
                    temp.legVec = temp.mPSIS-temp.mANK;
                    v = [0 0 1]; % Vertical line
%                     for i = 1:length(temp.u)
%                         TrialData(n).Results.thoraxObliq(i) = getAngVec(temp.u(i,:),v,2);
%                         TrialData(n).Results.legObliq(i) = getAngVec(temp.legVec(i,:),v,2);
%                     end
                    temp = [];
                end
                temp = [];
                % SWAY is calculated from the displacement of the Clavicle
                % markers, and referenced to the initial position.
                sway = sway - sway(1);
                COMsway = COMsway - COMsway(1);
                TrialData(n).Results.beamerSway = sway;
                TrialData(n).Results.clavY = clavY(start_idx:stop_idx);
                TrialData(n).Results.beamerCOMSway = COMsway;
                % Distance is calculated as the difference between the
                % forward positions of the clavicle marker at the beginning
                % and end of the analysis window.
                dist = clavY(stop_idx) - clavY(start_idx);
                if (TrialData(n).Info.Distance_Traveled == 144) 
                    TrialData(n).Results.totalDistance = 144*25.4-distOffset; % To be consistent with rest of code subracting out end of walking trial
                else
                    TrialData(n).Results.totalDistance = dist;
                end
                % Convert the distance to meters (from mm)
                TrialData(n).Results.totalDistance = TrialData(n).Results.totalDistance/1000;
                % We also store the average speed of the trial, and a
                % logical value indicating if the trial was completed.
                TrialData(n).Results.avgSpeed = TrialData(n).Results.totalDistance./time(end);
                TrialData(n).Results.completed = (TrialData(n).Info.Distance_Traveled == 144);
            end
            % RESULTS FOR ALL TRIALS
            TrialData(n).Results.time = time;
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