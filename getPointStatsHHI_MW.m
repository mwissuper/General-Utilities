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
    'meanPosPowerIntPtX','meanNegPowerIntPtX',... % Power at int. pt.
    'meanPosPowerIntPtY','meanNegPowerIntPtY',...
    'meanPosPowerIntPtZ','meanNegPowerIntPtZ',...
    'perPposFposX','perPnegFposX','perPposFnegX','perPnegFnegX',... % percentage of trial with certain combo's of sign in P and F
    'meanAbsPowerIntPtX','meanAbsPowerIntPtY','meanAbsPowerIntPtZ',... 
    'SDAbsPowerIntPtX','SDAbsPowerIntPtY','SDAbsPowerIntPtZ',...
    'meanAbsPowerPOBX','SDAbsPowerPOBX','meanPosPowerPOBX','meanNegPowerPOBX',... % Power as force at IP x velocity of torso - use this not for energy analysis but to understand haptic comm
    'POBperOppDevX','POBperAmpDevX','POBperOppRetX','POBperAmpRetX',... % percentage trial doing each type of action
    'meanPosFx','meanNegFx','SDPosFx','SDNegFx',... % Force metrics
    'meanPosFy','meanNegFy','SDPosFy','SDNegFy',...
    'meanPosFz','meanNegFz','SDPosFz','SDNegFz',...
    'meanFx','meanFy','meanFz',... % mean of force mag (abs force)
    'SDFx','SDFy','SDFz',... % SD of abs force 
    'rPowerFresX',... % correlation power and force model residual in x/ML dir
    'meanPosVx','meanNegVx','meanVx',... % Vel int pt metrics
    'meanPosVy','meanNegVy','meanVy',...
    'meanPosVz','meanNegVz','meanVz',...
    'meanVx_torso','meanVy_torso','meanVz_torso',... % Vel clav metrics
    'mx_torso','bx_torso','kx_torso','mx_lag_torso','bx_lag_torso','kx_lag_torso',...% regression coeff's of rFin marker acc, vel, pos to force in each dir
    'Rsqx_torso','Rsqx_lag_torso',...% Rsq of fit of regression above
    'Rsqx_IP','mx_IP','bx_IP','kx_IP',... % R^2 and reg coeff's IP state
    'lagFIPvTorsoX','xcorrFIPvTorsoX','lagvIPvTorsoX','xcorrvIPvTorsoX','xcorrFIPTorsoX','lagFIPTorsoX',...% xcorr of IP signals to POB Torso signals to see if hand interaction lags/leads balance changes
    'meanArmPOBX','SDarmPOBX'... % Magnitude of vector from 
    'StdSway','Dist','AvgSpeed',...% Other kinem metrics
    'Complete','Notes'}; 

% Old metrics
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
% Defaults are NaN for numeric entries, '' for string entries, and TRUE
% for logical entries
blank = {nan,nan,'',...
    nan,nan,...% power IP
    nan,nan,...
    nan,nan,...
    nan,nan,nan,nan,...% percentage of trial with certain combo's of sign in P (at IP) and F
    nan,nan,nan,...
    nan,nan,nan,...
    nan,nan,nan,nan,...% Power as force at IP x velocity of torso
    nan,nan,nan,nan,...% percentage of trial with certain combo's of sign in P (torso) and F
    nan,nan,nan,nan,... % force
    nan,nan,nan,nan,...
    nan,nan,nan,nan,...
    nan,nan,nan,...% abs force
    nan,nan,nan,...% net force
    nan,...% corr power and force model resid
    nan,nan,nan,...% vel IP
    nan,nan,nan,...
    nan,nan,nan,...
    nan,nan,nan,...% vel Torso
    nan,nan,nan,nan,nan,nan,...% regression Torso
    nan,nan,...% R^2 regression Torso
    nan,nan,nan,nan,... % R^2 and reg coeff's IP state
    nan,nan,nan,nan,nan,nan,...% xcorr of IP signals to POB Torso signals to see if hand interaction lags/leads balance changes
    nan,nan,...% POB arm "length" mean and SD in ML dir
    nan,nan,nan,...% other kinem's
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
            force = TrialData(n).Results.Forces; % This is with sign convention adjusted so compression is > 0 for Fz!
%             rFin = TrialData(n).Results.IntPt; % Use to get velocity of intPt
            % Old code arm len
%             [stats.PeakPosWork(n),stats.PeakNegWork(n)] = getPeaks(work);
%             [stats.PeakPosPower(n),stats.PeakNegPower(n)] = getPeaks(power);
            
            %% Interaction point power
            % Look at mean abs power and also mean pos and neg power. Also look
            % at proportion of pos power period when F > 0 and proportion of
            % neg power period when F < 0
            % Power   
            if ~isnan(TrialData(n).Results.IntPower) % bad trials
                for i = 1:3
                    clear temp indPos indNeg PosPower NegPower indPposFpos indPnegFpos
                    indPos = find(TrialData(n).Results.IntPower(:,i) > 0);
                    indNeg = find(TrialData(n).Results.IntPower(:,i) < 0);
                    PosPower = TrialData(n).Results.IntPower(indPos,i);
                    NegPower = TrialData(n).Results.IntPower(indNeg,i);
                    indPposFpos = find(TrialData(n).Results.Forces(indNeg,i) > 0);
                    indPnegFpos = find(TrialData(n).Results.Forces(indPos,i) > 0);
                    % Calculate these as a check
                    indPposFneg = find(TrialData(n).Results.Forces(indNeg,i) < 0);
                    indPnegFneg = find(TrialData(n).Results.Forces(indPos,i) < 0);
                    if i == 1
                        stats.meanAbsPowerIntPtX(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.SDAbsPowerIntPtX(n) = nanstd(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.meanPosPowerIntPtX(n) = nanmean(PosPower);
                        stats.meanNegPowerIntPtX(n) = nanmean(NegPower);
                        stats.perPposFposX(n) = length(indPposFpos)/length(TrialData(n).Results.IntPower(:,i));
                        stats.perPnegFposX(n) = length(indPnegFpos)/length(TrialData(n).Results.IntPower(:,i));
                        stats.perPposFnegX(n) = length(indPposFneg)/length(TrialData(n).Results.IntPower(:,i));
                        stats.perPnegFnegX(n) = length(indPnegFneg)/length(TrialData(n).Results.IntPower(:,i));
                    elseif i == 2
                        stats.meanAbsPowerIntPtY(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.SDAbsPowerIntPtY(n) = nanstd(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.meanPosPowerIntPtY(n) = nanmean(PosPower);
                        stats.meanNegPowerIntPtY(n) = nanmean(NegPower);
                    else
                        stats.meanAbsPowerIntPtZ(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.SDAbsPowerIntPtZ(n) = nanstd(abs(TrialData(n).Results.IntPower(:,i)));
                        stats.meanPosPowerIntPtZ(n) = nanmean(PosPower);
                        stats.meanNegPowerIntPtZ(n) = nanmean(NegPower);
                    end
                end

                %% POB power
                % Look at mean abs power and also mean pos and neg power
                % Power            
                clear temp indPos indNeg PosPOBpower NegPOBpower indPposFpos indPnegFpos indPposFneg indPnegFneg
                temp.POBpower = TrialData(n).Results.vTorso(:,1).*force(2:end,1);
                indPos = find(temp.POBpower > 0); % tension and move to right or compression and move to left (oppose dir of movement)
                indNeg = find(temp.POBpower < 0); % IP force amplifies direction of movement
                PosPOBpower = temp.POBpower(indPos);
                NegPOBpower = temp.POBpower(indNeg);

                stats.meanAbsPowerPOBX(n) = nanmean(abs(temp.POBpower));
                stats.SDAbsPowerPOBX(n) = nanstd(abs(temp.POBpower));
                stats.meanPosPowerPOBX(n) = nanmean(PosPOBpower);
                stats.meanNegPowerPOBX(n) = nanmean(NegPOBpower);
                
                % Percent time strategy
                indPposFpos = find(force(indNeg,1) > 0);
                indPnegFpos = find(force(indPos,1) > 0);
                indPposFneg = find(force(indNeg,1) < 0);
                indPnegFneg = find(force(indPos,1) < 0);
%                 stats.POBperPposFposX(n) = length(indPposFpos)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPnegFposX(n) = length(indPnegFpos)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPposFnegX(n) = length(indPposFneg)/length(TrialData(n).Results.vTorso(:,1));
%                 stats.POBperPnegFnegX(n) = length(indPnegFneg)/length(TrialData(n).Results.vTorso(:,1));
                
                %% POB x*x' and power
                % Look at relation of forces to person's movement,
                % reference of away or toward midline
                if strcmpi(stats.Type(n),'Assist Beam')
                    temp.midline = TrialData(n).Results.beamMidline;
                elseif strcmpi(stats.Type(n),'Assist Ground')
                    temp.midline = mean(TrialData(n).Results.torso(2:end,1));
                end
                temp.POBxv =  (TrialData(n).Results.torso(2:end,1)-temp.midline).*TrialData(n).Results.vTorso(:,1);
                indDev = find(temp.POBxv > 0); % deviation
                indindRet = find(temp.POBxv < 0); % return
                indOppDev = intersect(indDev,indPos); 
                indAmpDev = intersect(indDev,indNeg);
                indOppRet = intersect(indindRet,indPos); 
                indAmpRet = intersect(indindRet,indNeg);
                
                % Percent time strategy. All 4 should add up to 100%
                stats.POBperOppDevX(n) = length(indOppDev)/length(TrialData(n).Results.vTorso(:,1));
                stats.POBperAmpDevX(n) = length(indAmpDev)/length(TrialData(n).Results.vTorso(:,1));
                stats.POBperOppRetX(n) = length(indOppRet)/length(TrialData(n).Results.vTorso(:,1));
                stats.POBperAmpRetX(n) = length(indAmpRet)/length(TrialData(n).Results.vTorso(:,1));
                
                    
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

                %% Force - using lab/Vicon CS 1/31/20
                % Get mean (to compare to other studies) and SD of force for
                % pos, neg, net
                for i = 1:3
                    indPos = find(force(:,i) > 0);
                    indNeg = find(force(:,i) < 0);
                    if i == 1
                        stats.meanFx(n) = nanmean(abs(force(:,i)));
                        stats.SDFx(n) = nanstd(abs(force(:,i)));
                        stats.meanPosFx(n) = nanmean(force(indPos,i));
                        stats.SDPosFx(n) = nanstd(force(indPos,i));
                        stats.meanNegFx(n) = nanmean(force(indNeg,i));
                        stats.SDNegFx(n) = nanstd(force(indNeg,i));
                    elseif i == 2
                        stats.meanFy(n) = nanmean(abs(force(:,i)));
                        stats.SDFy(n) = nanstd(abs(force(:,i)));
                        stats.meanPosFy(n) = nanmean(force(indPos,i));
                        stats.SDPosFy(n) = nanstd(force(indPos,i));
                        stats.meanNegFy(n) = nanmean(force(indNeg,i));
                        stats.SDNegFy(n) = nanstd(force(indNeg,i));
                    else
                        stats.meanFz(n) = nanmean(abs(force(:,i)));
                        stats.SDFz(n) = nanstd(abs(force(:,i)));
                        stats.meanPosFz(n) = nanmean(force(indPos,i));
                        stats.SDPosFz(n) = nanstd(force(indPos,i));
                        stats.meanNegFz(n) = nanmean(force(indNeg,i));
                        stats.SDNegFz(n) = nanstd(force(indNeg,i));
                    end
                    clear indPos indNeg
                end

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

                %% Velocity metric to help interpret work data later
                % Get mean for pos, neg, abs
                for i = 1:3 
                    v = diff(TrialData(n).Results.IntPt(:,i)).*TrialData(n).Markers.samplerate;
                    indPos = find(v > 0);
                    indNeg = find(v < 0);
                    if i == 1
                        stats.meanVx(n) = nanmean(abs(v));
                        stats.meanPosVx(n) = nanmean(v(indPos));
                        stats.meanNegVx(n) = nanmean(v(indNeg));
                        stats.meanVx_torso(n) = nanmean(abs(TrialData(n).Results.vTorso(:,i)));
                    elseif i == 2
                        stats.meanVy(n) = nanmean(abs(v));
                        stats.meanPosVy(n) = nanmean(v(indPos));
                        stats.meanNegVy(n) = nanmean(v(indNeg));
                        stats.meanVy_torso(n) = nanmean(abs(TrialData(n).Results.vTorso(:,i)));
                    else
                        stats.meanVz(n) = nanmean(abs(v));
                        stats.meanPosVz(n) = nanmean(v(indPos));
                        stats.meanNegVz(n) = nanmean(v(indNeg));
                        stats.meanVz_torso(n) = nanmean(abs(TrialData(n).Results.vTorso(:,i)));
                    end
                    clear indPos indNeg v
                end

                %% Metrics for regression of POB state to force
                if TrialData(n).Results.px_torso < 0.05 % if fit is sig.
                    stats.mx_torso(n) = TrialData(n).Results.cx_torso(1);
                    stats.bx_torso(n) = TrialData(n).Results.cx_torso(2);
                    stats.kx_torso(n) = TrialData(n).Results.cx_torso(3);
                    stats.Rsqx_torso(n) = TrialData(n).Results.rsqx_torso;
                    % with lag
                    stats.mx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(1);
                    stats.bx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(2);
                    stats.kx_lag_torso(n) = TrialData(n).Results.cx_lag_torso(3);
                    stats.Rsqx_lag_torso(n) = TrialData(n).Results.rsqx_lag_torso;
                end
%                 if TrialData(n).Results.py_torso < 0.05 % if fit is sig.
%                     stats.my_torso(n) = TrialData(n).Results.cy_torso(1);
%                     stats.by_torso(n) = TrialData(n).Results.cy_torso(2);
%                     stats.ky_torso(n) = TrialData(n).Results.cy_torso(3);
%                     stats.Rsqy_torso(n) = TrialData(n).Results.rsqy_torso;
%                 end
%                 if TrialData(n).Results.pz_torso < 0.05 % if fit is sig.
%                     stats.mz_torso(n) = TrialData(n).Results.cz_torso(1);
%                     stats.bz_torso(n) = TrialData(n).Results.cz_torso(2);
%                     stats.kz_torso(n) = TrialData(n).Results.cz_torso(3);
%                     stats.Rsqz_torso(n) = TrialData(n).Results.rsqz_torso;
%                 end

                %% Metrics for regression of IP state to force
                if TrialData(n).Results.px_IP < 0.05 % if fit is sig.
                    stats.mx_IP(n) = TrialData(n).Results.cx_IP(1);
                    stats.bx_IP(n) = TrialData(n).Results.cx_IP(2);
                    stats.kx_IP(n) = TrialData(n).Results.cx_IP(3);
                    stats.Rsqx_IP(n) = TrialData(n).Results.rsqx_IP;
                end

                %% Correlation of power to regression residual

                % X/ML dir
                stats.rPowerFresX(n) = TrialData(n).Results.corrPowerFresX;

                %% Additional metrics that depend on force data

                %% Cross-correlation of IP signals to POB torso signals to see if interaction lags/leads balance (ML dir only)
                stats.xcorrFIPTorsoX(n) = TrialData(n).Results.xcorrFIPTorsoX;
                stats.lagFIPTorsoX(n) = TrialData(n).Results.lagFIPTorsoX;
                stats.xcorrFIPvTorsoX(n) = TrialData(n).Results.xcorrFIPvTorsoX;
                stats.lagFIPvTorsoX(n) = TrialData(n).Results.lagFIPvTorsoX;
                stats.xcorrvIPvTorsoX(n) = TrialData(n).Results.xcorrvIPvTorsoX;
                stats.lagvIPvTorsoX(n) = TrialData(n).Results.lagvIPvTorsoX;

                %% Correlation force to torso or COM pos
                stats.rFlatCOM(n) = TrialData(n).Results.rFlatCOM; 
                stats.rFsway(n) = TrialData(n).Results.rFsway;
                % Plots for all trials for a participant
                if plotoption > 0
                    plotind = plotind + 1;
                    colors = hsv(3);

                    if plotoption == 1 % Plot vector-based power with peaks
    %                     time = (0:(length(powerIntPtTot)-1))./TrialData(n).Markers.samplerate;
    %                     subplot(4,4,plotind),plot(time,powerIntPtTot),hold on, plot(indPkPosPowerTot,powerIntPtTot(indPkPosPowerTot),'o'),...
    %                         plot(indPkNegPowerTot,-powerIntPtTot(indPkNegPowerTot),'x')
    %                     legend('Power','Pos pks','Neg pks');
                    elseif plotoption == 2 % Plot indiv component power with peaks
                        time = (0:(length(powerIntPt(1,:))-1))./TrialData(n).Markers.samplerate;
                        subplot(4,5,plotind)
                            for i = 1:3
                                plot(powerIntPt(:,i),'color',colors(:,i)),hold on;
                            end
                            if plotind == 1
                                legend('x','y','z');
                            end
                            for i = 1:3 % Separate loop for legend's sake
                                temp = indPkPosPower{i};
                                if ~isnan(temp) % if there are non-nan values, i.e. if peaks exist
                                    plot(indPkPosPower{i},powerIntPt(i,temp),'o','color',colors(:,i)),...
                                end
                                temp = indPkNegPower{i};
                                if ~isnan(temp)
                                    plot(indPkNegPower{i},powerIntPt(i,temp),'x','color',colors(:,i))
                                end
                            end
                            hline(0,'k--');
                    elseif plotoption == 3 % Plot indiv components force with peaks
                        time = (0:(length(force(1,:))-1))./TrialData(n).Markers.samplerate;
                        subplot(4,5,plotind)
                            for i = 1:3
                                plot(force(:,i),'color',colors(:,i)),hold on;
                            end
                            if plotind == 1
                                legend('x','y','z');
                            end
                            for i = 1:3 % Separate loop for legend's sake
                                temp = indPkPosF{i};
                                if ~isnan(temp)
                                    plot(indPkPosF{i},force(i,temp),'o','color',colors(:,i)),...
                                end
                                temp = indPkNegF{i};
                                if ~isnan(temp)
                                    plot(indPkNegF{i},force(i,temp),'x','color',colors(:,i))
                                end
                            end
                            hline(0,'k--');
                    end
                    titlename = sprintf('%s %s',TrialData(n).Info.Trial,TrialData(n).Info.Condition);
                    xlabel('Frame'),ylabel('Force (N)'),title(titlename);
                end
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
