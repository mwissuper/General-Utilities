function stats = getPointStatsHHI_MW(TrialData,plotoption)
%% getPointStatsHHI: Point Statistics for the Human-Human Interaction Experiment
% 
%   STATS = getPointStatsHHI(TrialData) creates a table of point statistics
%   STATS from a structure array of experimental trial data TRIALDATA.
%   TrialData must contain a field INFO and a field Results. Ideally,
%   TrialData is the output of the function processHHI2017
%
%   STATS is a table containing point statistics for each trial in
%   TrialData. STATS also contains information about each trial, including
%   the subject number, trial number, trial type, and notes.

%   h is handle for figure
%   Luke Drnach
%   November 7, 2017

%   Refactored on December 5, 2018, for code clarity and reusability.%
%
%   Updated from HHI2017stats on December 5, 2018. All statistics are now 
%   stored in one table

%   TrialData is indexed by either early or late trial numbers so n =
%   1:numTrials in code below won't correspond to actual trial numbers

%   Calculate mean of peaks for power and total pos/neg work done for a
%   trial. MW 11/20/19.

%% Initialize the table with default values
% Write out the column names
varNames = {'Subject','TrialNumber','Type',...
    'meanVx','meanVz','meanVxz',... % RMS vel metrics
    'meanFx','meanFy','meanFz','meanFxz',... % mean of signed force
    'SDFx','SDFy','SDFz','SDFxz',... % SD of signed force
    'meanIPpowerX','meanIPpowerZ','meanIPpowerXZ','meanPOBpowerX','meanPOBpowerZ','meanPOBpowerXZ',... % mean of signed power at IP
    'SDIPpowerX','SDIPpowerZ','SDIPpowerXZ','SDPOBpowerX','SDPOBpowerZ','SDPOBpowerXZ',... % SD of signed power, similar to an RMS but not have sign issues 
    'meanTheta','SDTheta',... % Angle of force in XZ plane
    'mx_torso','bx_torso','kx_torso','Rsqx_torso',...% regression coeff's and fit of POB torso marker acc, vel, pos to force in lateral dir, with and without lag
    'mz_torso','bz_torso','kz_torso','Rsqz_torso',...% regression coeff's and R^2 of POB torso marker acc, vel, pos to force in vert dir
    'meanArmPOBX','SDarmPOBX'... % Magnitude of vector
    'StdSway','Dist','AvgSpeed',...% Other kinem metrics
    'meanLy','SDLy',... % Angular momentum mean and SD
    'meanTyGd','meanRMSTyGd','SDTyGd',... % mean and SD torque
    'meanTyFx','meanTyFz','SDTyFx','SDTyFz',...% torque components due to Fx and Fz
    'meanW','SDW','meanPang','SDPang',...% SD and signed mean angular vel and power
    'meanPangFx','meanPangFz','SDPangFx','SDPangFz',...% power components due to Fx and Fz
    'm_ang_torso','b_ang_torso','k_ang_torso','rsq_ang_torso',...% regression coeff's of POB torso marker acc, vel, pos to force in lateral dir, with and without lag
    'lagTyAngTorso','xcorrTyAngTorso','lagTyWTorso','xcorrTyWTorso',... % Xcorr Ty to torso ang state
    'Complete','Notes'}; 

%% Old metrics
    % 'mx_torso_st','bx_torso_st','kx_torso_st',... % Standardized coeff's from regression to normalized X
%     'mx_lag_torso','bx_lag_torso','kx_lag_torso','Rsqx_lag_torso',
%     'rPowerFresX',... % correlation power and force model residual in x/ML dir
%     'lagFIPvTorsoX','xcorrFIPvTorsoX','lagvIPvTorsoX','xcorrvIPvTorsoX','xcorrFIPTorsoX','lagFIPTorsoX',...% xcorr of IP signals to POB Torso signals to see if hand interaction lags/leads balance changes
%     'meanTorqueYTorso','SDTorqueYTorso',
%     'meanPosVx','meanNegVx','meanVx','SDVx',... % IP vel metrics
%     'meanPosVy','meanNegVy','meanVy',...
%     'meanPosVz','meanNegVz','meanVz','SDVz',...
%     'meanVx_torso','meanVy_torso','meanVz_torso',... % Vel clav metrics
%     'meanPosFx','meanNegFx','SDPosFx','SDNegFx',... 
%     'meanPosFy','meanNegFy','SDPosFy','SDNegFy',...
%     'meanPosFz','meanNegFz','SDPosFz','SDNegFz',...
%     'SDAbsPowerIntPtX','SDAbsPowerIntPtY','SDAbsPowerIntPtZ',...
%     'meanAbsPowerPOBX','SDAbsPowerPOBX','meanPosPowerPOBX','meanNegPowerPOBX',... % Power as force at IP x velocity of torso - use this not for energy analysis but to understand haptic comm
%     'POBperOppDevX','POBperAmpDevX','POBperOppRetX','POBperAmpRetX',... % percentage trial doing each type of action
%     'meanIPpowerPosF','meanIPpowerNegF',... % mean Pos and neg power for regular offset F drift (calc as check on IntPtPowerX)
%     'meanIPpowerPosFlo','meanIPpowerNegFlo',... % mean Pos and neg power for lo-boundary F drift
%     'meanIPpowerPosFhi','meanIPpowerNegFhi',... % mean Pos and neg power for hi-boundary F drift 
%     'perPposFposX','perPnegFposX','perPposFnegX','perPnegFnegX',... % percentage of trial with certain combo's of sign in P and F
%     'perVposFposX','perVnegFposX','perVposFnegX','perVnegFnegX',... % percentage of trial with certain combo's of sign in v and F
%     'perVposFloposX','perVnegFloposX','perVposFlonegX','perVnegFlonegX',... % percentage of trial with certain combo's of sign in v and lo-boundary drift F
%     'perVposFhiposX','perVnegFhiposX','perVposFhinegX','perVnegFhinegX',... % percentage of trial with certain combo's of sign in v and hi-boundary drift F
%     'meanPosPowerIntPtX','meanNegPowerIntPtX',... % Abs power at int. pt.
%     'meanPosPowerIntPtY','meanNegPowerIntPtY',...
%     'meanPosPowerIntPtZ','meanNegPowerIntPtZ',...
%     'POBperPposFposX','POBperPnegFposX','POBperPposFnegX','POBperPnegFnegX',... % percentage of trial with certain combo's of sign in P and F
%     'meanPosPowerPOBX','meanNegPowerPOBX',...
%     'meanPosPowerPOBY','meanNegPowerPOBY',...
%     'meanPosPowerPOBZ','meanNegPowerPOBZ',...
%     'meanAbsPowerPOBX','meanAbsPowerPOBY','meanAbsPowerPOBZ',...
%     'xcorrX','lagX','xcorrZ','lagZ',...% xcorr sway to force
%     'mx','bx','kx','my','by','ky','mz','bz','kz',...% regression coeff's of rFin marker acc, vel, pos to force in each dir
%     'Rsqx','Rsqy','Rsqz',...% Rsq of fit of regression above
%     'rTorsoFvert','rCOMFlat',
%     'PosWorkIntPtX','NegWorkIntPtX',...% Total work for trial (sum)
%     'PosWorkIntPtY','NegWorkIntPtY',...
%     'PosWorkIntPtZ','NegWorkIntPtZ',...
%     'meanPosWorkIntPtX','meanNegWorkIntPtX',...% Mean pos and neg work for trial per direction
%     'meanPosWorkIntPtY','meanNegWorkIntPtY',...
%     'meanPosWorkIntPtZ','meanNegWorkIntPtZ',...
%     'meanWorkIntPtX','meanWorkIntPtY','meanWorkIntPtZ',... % mean of net work (pos and neg combined)
%     'meanWorkMagIntPtX','meanWorkMagIntPtY','meanWorkMagIntPtZ',... % mean of abs work
%     'StdCOMSway',

%% Defaults are NaN for numeric entries, '' for string entries, and TRUE
% for logical entries
blank = {nan,nan,'',...
    nan,nan,nan,... % Vel metrics
    nan,nan,nan,nan,...% mean signed force
    nan,nan,nan,nan,...% SD signed force
    nan,nan,nan,nan,nan,nan,... % mean of signed power,
    nan,nan,nan,nan,nan,nan,... % SD of signed power
    nan,nan,... % Angle of force in XZ plane
    nan,nan,nan,nan,...% R^2 and reg coeff's POB in x
    nan,nan,nan,nan,... % R^2 and reg coeff's POB in z
    nan,nan,...% POB arm "length" mean and SD in ML dir
    nan,nan,nan,...% other kinem's
    nan,nan,... % Angular momentum mean and SD
    nan,nan,nan,...% mean and SD torque
    nan,nan,nan,nan,...% torque components due to Fx and Fz
    nan,nan,nan,nan,...% SD and mean signed angular vel and power
    nan,nan,nan,nan,...% power components due to Fx and Fz
    nan,nan,nan,nan,...% regression coeff's of POB torso marker acc, vel, pos to force in lateral dir, with and without lag
    nan,nan,nan,nan,...% Xcorr Ty to torso ang state
    true,' '};
% Filter out the STAT trials
Info = [TrialData.Info];
Trials = {Info.Trial};
statIdx = strncmp(Trials,'Stat',4);
TrialData(statIdx) = [];
% Copy the blank table for as many trials as there are
numTrials = length(TrialData);
blank = repmat(blank,numTrials,1);
% Create an 'empty' table
stats = cell2table(blank,'VariableNames',varNames);

% Butterworth 3rd order lowpass filter, cutoff 10Hz
[z, p, k] = butter(3,10/50);
[sos, g] = zp2sos(z, p, k);

plotind = 0;
%% Loop over the trials, pulling the point statistics
for n = 1:numTrials
    % Get the subject number 
    subject = TrialData(1).Info.Subject_1;
    idx = find(subject=='_') + 1;
    stats.Subject(n) = str2double(subject(idx:end));
    % Copy the trial information into the table    
    stats.TrialNumber(n) = str2double(TrialData(n).Info.Trial(end-1:end));
    stats.Type{n} = TrialData(n).Info.Condition;
    if ~strcmpi(stats.Type(n),'Assist Solo') & ~strcmpi(stats.Type(n),'Assist Solos') % No POB
        %% There are no values to copy for Assist Solo type trials. For all
        % trials, we add in the standard deviation of the POB sway, the
        % distance traveled, the average speed, and the completion status.
        stats.StdSway(n) = std(TrialData(n).Results.beamerSway);
        stats.StdCOMSway(n) = std(TrialData(n).Results.beamerCOMSway);
%         stats.StdPelvicObliq(n) = std(TrialData(n).Results.pelvicObliq);
%         stats.StdThoraxObliq(n) = std(TrialData(n).Results.thoraxObliq);
%         stats.StdLegObliq(n) = std(TrialData(n).Results.legObliq);
%         % Calculate correlation two obliquity angles to see if in or out of
%         % phase. Use rho value only if statistically sig.
%         [rho stats.pPelvicThoraxObliq(n)] = corr(TrialData(n).Results.pelvicObliq,TrialData(n).Results.thoraxObliq');
%         if stats.pPelvicThoraxObliq(n) < 0.05
%             stats.rPelvicThoraxObliq(n) = rho;
%         else
%             stats.rPelvicThoraxObliq(n) = nan;
%         end
        % Calculate correlation two obliquity angles to see if in or out of
        % phase. Use rho value only if statistically sig.
%         [rho stats.pLegThoraxObliq(n)] = corr(TrialData(n).Results.legObliq',TrialData(n).Results.thoraxObliq');
%         if stats.pLegThoraxObliq(n) < 0.05
%             stats.rLegThoraxObliq(n) = rho;
%         else
%             stats.rLegThoraxObliq(n) = nan;
%         end
        stats.Dist(n) = TrialData(n).Results.totalDistance;
        stats.AvgSpeed(n) = TrialData(n).Results.avgSpeed;
        stats.Complete(n) = TrialData(n).Results.completed;    
        % Torso angular speed (mean RMS)
        stats.meanW(n) = nanmean(sqrt(TrialData(n).Results.wTorso.^2));
        stats.SDW(n) = nanstd(sqrt(TrialData(n).Results.wTorso.^2));
        
        stats.meanLy(n) = nanmean(sqrt(TrialData(n).Results.Ly.^2));
        stats.SDLy(n) = nanstd(sqrt(TrialData(n).Results.Ly.^2));
        if any(strcmpi(stats.Type(n),{'Assist Ground','Assist Beam'}))
            %% For Assisted trials, we must calculate peak work, power, and
            % force, and add it to the stats table
            %Unpack values from the Results field
%             power = TrialData(n).Results.AssistPower; % calculated per dir from multiply force by velocity finger
%             work = TrialData(n).Results.AssistCumWork; % calculated per dir from integrate power
%             powerIntPtTot = TrialData(n).Results.IntPower_tot;
%             workIntPtTot = TrialData(n).Results.IntCumWork_tot;
%             powerIntPt = TrialData(n).Results.IntPower;
%             workIntPt = TrialData(n).Results.IntCumWork; % work calculated along each dir separately
            force = TrialData(n).Results.Force; % This is with sign convention adjusted to be force on POB
%             rFin = TrialData(n).Results.IntPt; % Use to get velocity of intPt
            % Old code arm len
%             [stats.PeakPosWork(n),stats.PeakNegWork(n)] = getPeaks(work);
%             [stats.PeakPosPower(n),stats.PeakNegPower(n)] = getPeaks(power);
            
            %% Interaction point power
            
            % Look at force and force lo/hi bounds with drift and combo
            % with velocity
                        
            % Look at mean abs power and also mean pos and neg power. Also look
            % at proportion of pos power period when F > 0 and proportion of
            % neg power period when F < 0
 
            if ~isnan(TrialData(n).Results.IntPower) % good trials                                           
                %% Force - using lab/Vicon CS 1/31/20
                % Get signed mean (to compare to other studies) and SD of force for

                for i = 1:3
                    indPos = find(force(:,i) > 0);
                    indNeg = find(force(:,i) < 0);
                    if i == 1
                        stats.meanFx(n) = nanmean(force(:,i));
                        stats.SDFx(n) = nanstd(force(:,i));
%                         stats.meanPosFx(n) = nanmean(force(indPos,i));
%                         stats.SDPosFx(n) = nanstd(force(indPos,i));
%                         stats.meanNegFx(n) = nanmean(force(indNeg,i));
%                         stats.SDNegFx(n) = nanstd(force(indNeg,i));
                    elseif i == 2
                        stats.meanFy(n) = nanmean(force(:,i));
                        stats.SDFy(n) = nanstd(abs(force(:,i)));
%                         stats.meanPosFy(n) = nanmean(force(indPos,i));
%                         stats.SDPosFy(n) = nanstd(force(indPos,i));
%                         stats.meanNegFy(n) = nanmean(force(indNeg,i));
%                         stats.SDNegFy(n) = nanstd(force(indNeg,i));
                    else
                        stats.meanFz(n) = nanmean(force(:,i));
                        stats.SDFz(n) = nanstd(abs(force(:,i)));
%                         stats.meanPosFz(n) = nanmean(force(indPos,i));
%                         stats.SDPosFz(n) = nanstd(force(indPos,i));
%                         stats.meanNegFz(n) = nanmean(force(indNeg,i));
%                         stats.SDNegFz(n) = nanstd(force(indNeg,i));
                    end
                    clear indPos indNeg
                end
                
                % 2D vector force
                clear Fmag theta
                for i = 1:length(force(:,1))
                    Fmag(i) = norm(force(i,[1 3]));
                    theta(i) = atan2(force(i,3),force(i,1));
                end
                stats.meanFxz(n) = nanmean(Fmag);
                stats.SDFxz(n) = nanstd(Fmag);
                stats.meanTheta(n) = nanmean(theta);
                stats.SDTheta(n) = nanstd(theta);
    %             stats.PeakNegFx(n) = mean(PkNegFx);
                % X dir
    %             % Seems arbitrary to pick what counts as a peak and to divide
    %             % epochs
    %             % X dir
    %             [PkPosFx,indPkPosF{1},PkNegFx,indPkNegF{1}] = getPosNegPeaks(force(1,:));
    %             stats.PeakPosFx(n) = mean(PkPosFx);
    %             stats.PeakNegFx(n) = mean(PkNegFx);
    %             % X dir
    %             [PkPosFy,indPkPosF{2},PkNegFy,indPkNegF{2}] = getPosNegPeaks(force(2,:));
    %             stats.PeakPosFy(n) = mean(PkPosFy);
    %             stats.PeakNegFy(n) = mean(PkNegFy);
    %             % Z dir
    %             [PkZCompression,indPkPosF{3},PkZTension,indPkNegF{3}] = getPosNegPeaks(force(3,:));
    %             stats.PeakZCompression(n) = mean(PkZCompression);
    %             stats.PeakZTension(n) = mean(PkZTension);
    %             % Old code just looked at one max value for trial instead of
    %             mean of extrema/peaks
    %             [stats.PeakZTension(n),stats.PeakZCompression(n)] = getPeaks(force(3,:));
    %             [stats.PeakPosFy(n), stats.PeakNegFy(n)] = getPeaks(force(2,:));
    %             [stats.PeakPosFx(n),stats.PeakNegFx(n)] = getPeaks(force(1,:));
    
                %% IP Velocity metric to help interpret other metrics
                % Calculate mean of RMS lateral vel to compare to Dagmar's
                % work
                for i = 1:3 
                    v = diff(TrialData(n).Results.IntPt(:,i)).*TrialData(n).Markers.samplerate;
                    indPos = find(v > 0);
                    indNeg = find(v < 0);
                    if i == 1
%                         stats.SDVx(n) = nanstd(v);
%                         stats.meanVx(n) = nanmean(abs(v));
%                         stats.SDVx(n) = nanstd(abs(v));
                        stats.meanVx(n) = nanmean(sqrt(v.^2));
%                         stats.SDVx(n) = nanstd(sqrt(v.^2));
%                         stats.meanPosVx(n) = nanmean(v(indPos));
%                         stats.meanNegVx(n) = nanmean(v(indNeg));
%                         stats.meanVx_torso(n) = nanmean(sqrt(TrialData(n).Results.vTorso(:,i).^2));
                    elseif i == 2
%                         stats.meanVy(n) = nanmean(abs(v));
%                         stats.meanPosVy(n) = nanmean(v(indPos));
%                         stats.meanNegVy(n) = nanmean(v(indNeg));
%                         stats.meanVy_torso(n) = nanmean(sqrt(TrialData(n).Results.vTorso(:,i).^2));
                    else
%                         stats.SDVz(n) = nanstd(v);
                        stats.meanVz(n) = nanmean(sqrt(v.^2));
%                         stats.SDVz(n) = nanstd(sqrt(v.^2));
%                         stats.meanPosVz(n) = nanmean(v(indPos));
%                         stats.meanNegVz(n) = nanmean(v(indNeg));
%                         stats.meanVz_torso(n) = nanmean(sqrt(TrialData(n).Results.vTorso(:,i).^2));
                    end
                    clear indPos indNeg v
                end
                
                % 2D vector of IP
                clear vmag
                for i = 1:length(TrialData(n).Results.IntPtVel(:,1))
                    vmag(i) = norm(TrialData(n).Results.IntPtVel(i,[1 3]));
                end
                stats.meanVxz(n) = nanmean(vmag);
%                 stats.SDVxz(n) = nanstd(vmag);                

                %% Linear power. Mostly care SD of signed power. Don't look pos and neg portions as unsure force sensor drift
                for i = 1:3
                    clear temp indPos indNeg PosPower NegPower indPposFpos indPnegFpos
%                     indPos = find(TrialData(n).Results.IntPower(:,i) > 0);
%                     indNeg = find(TrialData(n).Results.IntPower(:,i) < 0);
%                     PosPower = TrialData(n).Results.IntPower(indPos,i);
%                     NegPower = TrialData(n).Results.IntPower(indNeg,i);
%                     indPposFpos = find(TrialData(n).Results.Force(indNeg,i) > 0);
%                     indPnegFpos = find(TrialData(n).Results.Force(indPos,i) > 0);
%                     % Calculate these as a check
%                     indPposFneg = find(TrialData(n).Results.Force(indNeg,i) < 0);
%                     indPnegFneg = find(TrialData(n).Results.Force(indPos,i) < 0);
                    if i == 1
                        stats.meanIPpowerX(n) = nanmean(TrialData(n).Results.IntPower(:,i));
                        stats.SDIPpowerX(n) = nanstd(TrialData(n).Results.IntPower(:,i));
                        stats.meanPOBpowerX(n) = nanmean(TrialData(n).Results.POBPower(:,i));
                        stats.SDPOBpowerX(n) = nanstd(TrialData(n).Results.POBPower(:,i));
%                         stats.meanPosPowerIntPtX(n) = nanmean(PosPower);
%                         stats.meanNegPowerIntPtX(n) = nanmean(NegPower);
%                         stats.perPposFposX(n) = length(indPposFpos)/length(TrialData(n).Results.IntPower(:,i));
%                         stats.perPnegFposX(n) = length(indPnegFpos)/length(TrialData(n).Results.IntPower(:,i));
%                         stats.perPposFnegX(n) = length(indPposFneg)/length(TrialData(n).Results.IntPower(:,i));
%                         stats.perPnegFnegX(n) = length(indPnegFneg)/length(TrialData(n).Results.IntPower(:,i));
                    elseif i == 2
%                         stats.meanIPpowerY(n) = nanmean(TrialData(n).Results.IntPower(:,i));
%                         stats.SDIPpowerY(n) = nanstd(TrialData(n).Results.IntPower(:,i));
%                         stats.meanPosPowerIntPtY(n) = nanmean(PosPower);
%                         stats.meanNegPowerIntPtY(n) = nanmean(NegPower);
                    else
                        stats.meanIPpowerZ(n) = nanmean(TrialData(n).Results.IntPower(:,i));
                        stats.SDIPpowerZ(n) = nanstd(TrialData(n).Results.IntPower(:,i));
                        stats.meanPOBpowerZ(n) = nanmean(TrialData(n).Results.POBPower(:,i));
                        stats.SDPOBpowerZ(n) = nanstd(TrialData(n).Results.POBPower(:,i));
%                         stats.meanPosPowerIntPtZ(n) = nanmean(PosPower);
%                         stats.meanNegPowerIntPtZ(n) = nanmean(NegPower);
                    end
                end
                
                % 2D power
                stats.meanIPpowerXZ(n) = nanmean(TrialData(n).Results.IntPowerXZ);
                stats.SDIPpowerXZ(n) = nanstd(TrialData(n).Results.IntPowerXZ);
                stats.meanPOBpowerXZ(n) = nanmean(TrialData(n).Results.POBPowerXZ);
                stats.SDPOBpowerXZ(n) = nanstd(TrialData(n).Results.POBPowerXZ);
                
                %% Torque metrics 
                stats.meanTyGd(n) = nanmean(TrialData(n).Results.TyGd); % Signed mean
                stats.meanRMSTyGd(n) = nanmean(sqrt(TrialData(n).Results.TyGd.^2)); % RMS mean to use to compare other studies
                stats.SDTyGd(n) = nanstd(TrialData(n).Results.TyGd); 
%                 stats.meanTorqueYTorso(n) = nanmean(sqrt(TrialData(n).Results.TyTorso.^2));
%                 stats.SDTorqueYTorso(n) = nanstd(sqrt(TrialData(n).Results.TyTorso.^2));
                stats.meanTyFx(n) = nanmean(TrialData(n).Results.TyFx); % Signed mean
                stats.SDTyFx(n) = nanstd(TrialData(n).Results.TyFx); 
                stats.meanTyFz(n) = nanmean(TrialData(n).Results.TyFz); % Signed mean
                stats.SDTyFz(n) = nanstd(TrialData(n).Results.TyFz); 
                                
                %% Angular power
                stats.meanPang(n) = nanmean(TrialData(n).Results.angP); % signed
                stats.SDPang(n) = nanstd(TrialData(n).Results.angP);
                stats.meanPangRMS(n) = nanmean(sqrt(TrialData(n).Results.angP.^2)); % signed
                stats.SDPangRMS(n) = nanstd(sqrt(TrialData(n).Results.angP.^2)); % signed
                stats.meanPangFx(n) = nanmean(TrialData(n).Results.angPFx); % signed
                stats.SDPangFx(n) = nanstd(TrialData(n).Results.angPFx);
                stats.meanPangFz(n) = nanmean(TrialData(n).Results.angPFz); % signed
                stats.SDPangFz(n) = nanstd(TrialData(n).Results.angPFz);
%                 stats.meanPowerAngRMS(n) = nanmean(sqrt(TrialData(n).Results.angP.^2));
%                 clear indPos indNeg
%                 indPos = find(TrialData(n).Results.angP > 0); 
%                 indNeg = find(TrialData(n).Results.angP < 0); 
%                 stats.meanPangPos(n) = nanmean(TrialData(n).Results.angP(indPos));
%                 stats.meanPangNeg(n) = nanmean(TrialData(n).Results.angP(indNeg));
                

                %% Correlation metrics (angular)
                stats.lagTyAngTorso(n) = nanmean(TrialData(n).Results.lagTyAngTorso);
                stats.xcorrTyAngTorso(n) = nanmean(TrialData(n).Results.xcorrTyAngTorso);
                stats.lagTyWTorso(n) = nanmean(TrialData(n).Results.lagTyWTorso);
                stats.xcorrTyWTorso(n) = nanmean(TrialData(n).Results.xcorrTyWTorso);
                % Don't trust alphaTorso

                %% Metrics for regression of POB state to force
                
                % ML/x dir
                if TrialData(n).Results.px_torso < 0.05 % if fit is sig.
                    stats.mx_torso(n) = TrialData(n).Results.cx_torso(1);
                    stats.bx_torso(n) = TrialData(n).Results.cx_torso(2);
                    stats.kx_torso(n) = TrialData(n).Results.cx_torso(3);
                    stats.Rsqx_torso(n) = TrialData(n).Results.rsqx_torso;
                    % Standardized version of above (R^2) should be the same
                    stats.mx_torso_st(n) = TrialData(n).Results.cx_torso_st(1);
                    stats.bx_torso_st(n) = TrialData(n).Results.cx_torso_st(2);
                    stats.kx_torso_st(n) = TrialData(n).Results.cx_torso_st(3);
                    stats.Rsqx_torso_st(n) = TrialData(n).Results.rsqx_torso_st;
%                     % with lag
%                     stats.mx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(1);
%                     stats.bx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(2);
%                     stats.kx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(3);
%                     stats.Rsqx_lag_torso(n) = TrialData(n).Results.rsqx_lag_torso;
                end
%                 if TrialData(n).Results.py_torso < 0.05 % if fit is sig.
%                     stats.my_torso(n) = TrialData(n).Results.cy_torso(1);
%                     stats.by_torso(n) = TrialData(n).Results.cy_torso(2);
%                     stats.ky_torso(n) = TrialData(n).Results.cy_torso(3);
%                     stats.Rsqy_torso(n) = TrialData(n).Results.rsqy_torso;
%                 end

                % Vert/z dir
                if TrialData(n).Results.pz_torso < 0.05 % if fit is sig.
                    stats.mz_torso(n) = TrialData(n).Results.cz_torso(1);
                    stats.bz_torso(n) = TrialData(n).Results.cz_torso(2);
                    stats.kz_torso(n) = TrialData(n).Results.cz_torso(3);
                    stats.Rsqz_torso(n) = TrialData(n).Results.rsqz_torso;
                end
                
                % Angular
                if TrialData(n).Results.p_ang_torso < 0.05 % if fit is sig.
                    stats.m_ang_torso(n) = TrialData(n).Results.c_ang_torso(1);
                    stats.b_ang_torso(n) = TrialData(n).Results.c_ang_torso(2);
                    stats.k_ang_torso(n) = TrialData(n).Results.c_ang_torso(3);
                    stats.Rsq_ang_torso(n) = TrialData(n).Results.rsq_ang_torso;
                end

                %% Correlation of power to regression residual

                % X/ML dir
%                 stats.rPowerFresX(n) = TrialData(n).Results.corrPowerFresX;

                %% Cross-correlation of IP signals to POB torso signals to see if interaction lags/leads balance (ML dir only)
%                 stats.xcorrFIPTorsoX(n) = TrialData(n).Results.xcorrFIPTorsoX;
%                 stats.lagFIPTorsoX(n) = TrialData(n).Results.lagFIPTorsoX;
%                 stats.xcorrFIPvTorsoX(n) = TrialData(n).Results.xcorrFIPvTorsoX;
%                 stats.lagFIPvTorsoX(n) = TrialData(n).Results.lagFIPvTorsoX;
%                 stats.xcorrvIPvTorsoX(n) = TrialData(n).Results.xcorrvIPvTorsoX;
%                 stats.lagvIPvTorsoX(n) = TrialData(n).Results.lagvIPvTorsoX;

                %% IP force and v IP combo's (x dir only)
%                 clear vx 
%                 vx = diff(filtfilt(sos,g,TrialData(n).Results.IntPt(:,1))).*TrialData(n).Markers.samplerate;
%                 [stats.perVposFposX(n),stats.perVnegFposX(n),stats.perVposFnegX(n),stats.perVnegFnegX(n),stats.meanIPpowerPosF(n),stats.meanIPpowerNegF(n)] = signProdPer(vx,TrialData(n).Results.Force(2:end,1)); % meanABpos should be same as mean of intPt_power_pos
%                 
%                 % Do same as above for lo boundary F with drift
%                 [stats.perVposFloposX(n),stats.perVnegFloposX(n),stats.perVposFlonegX(n),stats.perVnegFlonegX(n),stats.meanIPpowerPosFlo(n),stats.meanIPpowerNegFlo(n)] = signProdPer(vx,TrialData(n).Results.Force(2:end,1)-1.5); % meanABpos should be same as mean of intPt_power_pos
%                 
%                 % Do same as above for hi boundary F with drift
%                 [stats.perVposFhiposX(n),stats.perVnegFhiposX(n),stats.perVposFhinegX(n),stats.perVnegFhinegX(n),stats.meanIPpowerPosFhi(n),stats.meanIPpowerNegFhi(n)] = signProdPer(vx,TrialData(n).Results.Force(2:end,1)+1.5); % meanABpos should be same as mean of intPt_power_pos

%                 %% POB power
%                 % Look at mean abs power and also mean pos and neg power
%                 % Power            
%                 clear temp indPos indNeg PosPOBpower NegPOBpower indPposFpos indPnegFpos indPposFneg indPnegFneg
%                 temp.POBpower = TrialData(n).Results.vTorso(:,1).*force(2:end,1);
%                 indPos = find(temp.POBpower > 0); % tension and move to right or compression and move to left (oppose dir of movement)
%                 indNeg = find(temp.POBpower < 0); % IP force amplifies direction of movement
%                 PosPOBpower = temp.POBpower(indPos);
%                 NegPOBpower = temp.POBpower(indNeg);
% 
%                 stats.meanAbsPowerPOBX(n) = nanmean(abs(temp.POBpower));
%                 stats.SDAbsPowerPOBX(n) = nanstd(abs(temp.POBpower));
%                 stats.meanPosPowerPOBX(n) = nanmean(PosPOBpower);
%                 stats.meanNegPowerPOBX(n) = nanmean(NegPOBpower);
%                 
%                 % Percent time strategy
%                 indPposFpos = find(force(indNeg,1) > 0);
%                 indPnegFpos = find(force(indPos,1) > 0);
%                 indPposFneg = find(force(indNeg,1) < 0);
%                 indPnegFneg = find(force(indPos,1) < 0);
%                 stats.POBperPposFposX(n) = length(indPposFpos)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPnegFposX(n) = length(indPnegFpos)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPposFnegX(n) = length(indPposFneg)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPnegFnegX(n) = length(indPnegFneg)/length(TrialData(n).Results.vTorso(:,1));
                
%                 %% POB x*x' and power
%                 % Look at relation of forces to person's movement,
%                 % reference of away or toward midline
%                 if strcmpi(stats.Type(n),'Assist Beam')
%                     temp.midline = TrialData(n).Results.midline;
%                 elseif strcmpi(stats.Type(n),'Assist Ground')
%                     temp.midline = mean(TrialData(n).Results.torso(2:end,1));
%                 end
%                 temp.POBxv = (TrialData(n).Results.torso(2:end,1)-temp.midline).*TrialData(n).Results.vTorso(:,1);
%                 indDev = find(temp.POBxv > 0); % deviation
%                 indindRet = find(temp.POBxv < 0); % return
%                 indOppDev = intersect(indDev,indPos); 
%                 indAmpDev = intersect(indDev,indNeg);
%                 indOppRet = intersect(indindRet,indPos); 
%                 indAmpRet = intersect(indindRet,indNeg);
%                 
%                 % Percent time strategy. All 4 should add up to 100%
%                 stats.POBperOppDevX(n) = length(indOppDev)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperAmpDevX(n) = length(indAmpDev)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperOppRetX(n) = length(indOppRet)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperAmpRetX(n) = length(indAmpRet)/length(TrialData(n).Results.vTorso(:,1));
                                    
                %% Mean Pos and Neg Work done per trial
        % %             % Dot product vectors
        % %             [stats.PosWorkIntPtTot(n),stats.NegWorkIntPtTot(n),stats.meanPosWorkIntPtTot(n),stats.meanNegWorkIntPtTot(n)] = getPosNegWork(powerIntPtTot,TrialData(n).Markers.samplerate);
        %             % Each direction separately
        %             for i = 1:3
        %                [tempPos,tempNeg,tempPosMean,tempNegMean] = getPosNegWork(powerIntPt(:,i),TrialData(n).Markers.samplerate); % This calculates one number for work using periods of time where power was pos or neg and integrating
        %                 if i == 1
        %                     stats.PosWorkIntPtX(n) = tempPos; % Total work done for whole trial (cumsum evaluated at end of trial)
        %                     stats.NegWorkIntPtX(n) = tempNeg;
        %                     stats.meanPosWorkIntPtX(n) = tempPosMean; % Mean work done during trial
        %                     stats.meanNegWorkIntPtX(n) = tempNegMean;
        %                 elseif i == 2
        %                     stats.PosWorkIntPtY(n) = tempPos;
        %                     stats.NegWorkIntPtY(n) = tempNeg;
        %                     stats.meanPosWorkIntPtY(n) = tempPosMean;
        %                     stats.meanNegWorkIntPtY(n) = tempNegMean;
        %                 else
        %                     stats.PosWorkIntPtZ(n) = tempPos;
        %                     stats.NegWorkIntPtZ(n) = tempNeg;
        %                     stats.meanPosWorkIntPtZ(n) = tempPosMean;
        %                     stats.meanNegWorkIntPtZ(n) = tempNegMean;
        %                 end
        %             end
        %             
        %             %% Mean of net (pos + neg) and abs work vs. time per trial
        %             % Each direction separately
        %             for i = 1:3
        %                 if i == 1
        %                     stats.meanWorkIntPtX(n) = nanmean(workIntPt(:,i));
        %                     stats.meanWorkMagIntPtX(n) = nanmean(abs(workIntPt(:,i)));
        %                 elseif i == 2
        %                     stats.meanWorkIntPtY(n) = nanmean(workIntPt(:,i));
        %                     stats.meanWorkMagIntPtY(n) = nanmean(abs(workIntPt(:,i)));
        %                 else
        %                     stats.meanWorkIntPtZ(n) = nanmean(workIntPt(:,i));
        %                     stats.meanWorkMagIntPtZ(n) = nanmean(abs(workIntPt(:,i)));
        %                 end
        %             end

            end
            %% Distance torso to IP in ML dir
            stats.meanArmPOBX(n) = nanmean(TrialData(n).Results.armPOBX);
            stats.SDarmPOBX(n) = nanstd(TrialData(n).Results.armPOBX);
        end
        % Copy Notes to the table (if they exist)
        if isfield(TrialData(n).Info,'Notes')
            stats.Notes{n} = char(TrialData(n).Info.Notes);
        end
    end
end
% if nargout == 0
%     % If there are no output arguments, make plots of the statistics
%     plotHHIStats(stats);
% end
end

function [pmax,pmin] = getPeaks(data)
% GETPEAKS: Helper function calculating the minimum (below zero) and
% maximum (above zero) values.
pmax = max(data(data>0));
if isempty(pmax)
    pmax = NaN;
end
pmin = min(data(data<0));
if isempty(pmin)
    pmin = NaN;
end
end
