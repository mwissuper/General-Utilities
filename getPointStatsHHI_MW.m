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
    'meanAbsPowerIntPtX','meanAbsPowerIntPtY','meanAbsPowerIntPtZ',... 
    'meanPosFx','meanNegFx','SDPosFx','SDNegFx',... % Force metrics
    'meanPosFy','meanNegFy','SDPosFy','SDNegFy',...
    'meanPosFz','meanNegFz','SDPosFz','SDNegFz',...
    'meanFx','meanFy','meanFz',... % mean of force mag (abs force)
    'SDFx','SDFy','SDFz',... % SD of net force (pos and neg combined)
    'rPowerFresX',... % correlation power and force model residual in x/ML dir
    'meanPosVx','meanNegVx','meanVx',... % Vel int pt metrics
    'meanPosVy','meanNegVy','meanVy',...
    'meanPosVz','meanNegVz','meanVz',...
    'meanVx_clav','meanVy_clav','meanVz_clav',... % Vel clav metrics
    'mx_clav','bx_clav','kx_clav','my_clav','by_clav','ky_clav','mz_clav','bz_clav','kz_clav',...% regression coeff's of rFin marker acc, vel, pos to force in each dir
    'Rsqx_clav','Rsqy_clav','Rsqz_clav',...% Rsq of fit of regression above
    'StdSway','Dist','AvgSpeed',...% Other kinem metrics
    'Complete','Notes'}; 

% Old metrics
%     'meanPosPowerPOBX','meanNegPowerPOBX',...
%     'meanPosPowerPOBY','meanNegPowerPOBY',...
%     'meanPosPowerPOBZ','meanNegPowerPOBZ',...
%     'meanAbsPowerPOBX','meanAbsPowerPOBY','meanAbsPowerPOBZ',...
%     'xcorrX','lagX','xcorrZ','lagZ',...% xcorr sway to force
%     'mx','bx','kx','my','by','ky','mz','bz','kz',...% regression coeff's of rFin marker acc, vel, pos to force in each dir
%     'Rsqx','Rsqy','Rsqz',...% Rsq of fit of regression above
%     'rClavFvert','rCOMFlat',
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
    nan,nan,nan,...
    nan,nan,nan,nan,... % force
    nan,nan,nan,nan,...
    nan,nan,nan,nan,...
    nan,nan,nan,...% abs force
    nan,nan,nan,...% net force
    nan,...% corr power and force model resid
    nan,nan,nan,...% vel IP
    nan,nan,nan,...
    nan,nan,nan,...
    nan,nan,nan,...% vel clav
    nan,nan,nan,nan,nan,nan,nan,nan,nan,...% regression clav
    nan,nan,nan,...% R^2 regression clav
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
%         stats.StdCOMSway(n) = std(TrialData(n).Results.beamerCOMSway);
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
            rFin = TrialData(n).Results.IntPt; % Use to get velocity of intPt
            % Old code arm len
%             [stats.PeakPosWork(n),stats.PeakNegWork(n)] = getPeaks(work);
%             [stats.PeakPosPower(n),stats.PeakNegPower(n)] = getPeaks(power);
            
        %% Interaction point power
        % Look at mean abs power and also mean pos and neg power
        % Power            
        for i = 1:3
            clear temp indPos indNeg PosPower NegPower
            indPos = find(TrialData(n).Results.IntPower(:,i) > 0);
            indNeg = find(TrialData(n).Results.IntPower(:,i) < 0);
            PosPower = TrialData(n).Results.IntPower(indPos,i);
            NegPower = TrialData(n).Results.IntPower(indNeg,i);
            if i == 1
                stats.meanAbsPowerIntPtX(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                stats.meanPosPowerIntPtX(n) = nanmean(PosPower);
                stats.meanNegPowerIntPtX(n) = nanmean(NegPower);
            elseif i == 2
                stats.meanAbsPowerIntPtY(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                stats.meanPosPowerIntPtY(n) = nanmean(PosPower);
                stats.meanNegPowerIntPtY(n) = nanmean(NegPower);
            else
                stats.meanAbsPowerIntPtZ(n) = nanmean(abs(TrialData(n).Results.IntPower(:,i)));
                stats.meanPosPowerIntPtZ(n) = nanmean(PosPower);
                stats.meanNegPowerIntPtZ(n) = nanmean(NegPower);
            end
        end
        
        %% POB power
        % Look at mean abs power and also mean pos and neg power
        % Power            
%         for i = 1:3
%             clear temp indPos indNeg PosPOBpower NegPOBpower
%             indPos = find(TrialData(n).Results.POBpower(:,i) > 0);
%             indNeg = find(TrialData(n).Results.POBpower(:,i) < 0);
%             PosPOBpower = TrialData(n).Results.POBpower(indPos,i);
%             NegPOBpower = TrialData(n).Results.POBpower(indNeg,i);
%             if i == 1
%                 stats.meanAbsPowerPOBX(n) = nanmean(abs(TrialData(n).Results.POBpower(:,i)));
%                 stats.meanPosPowerPOBX(n) = nanmean(PosPOBpower);
%                 stats.meanNegPowerPOBX(n) = nanmean(NegPOBpower);
%             elseif i == 2
%                 stats.meanAbsPowerPOBY(n) = nanmean(abs(TrialData(n).Results.POBpower(:,i)));
%                 stats.meanPosPowerPOBY(n) = nanmean(PosPOBpower);
%                 stats.meanNegPowerPOBY(n) = nanmean(NegPOBpower);
%             else
%                 stats.meanAbsPowerPOBZ(n) = nanmean(abs(TrialData(n).Results.POBpower(:,i)));
%                 stats.meanPosPowerPOBZ(n) = nanmean(PosPOBpower);
%                 stats.meanNegPowerPOBZ(n) = nanmean(NegPOBpower);
%             end
%         end
%             
%             %% Mean Pos and Neg Work done per trial
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
                    stats.meanVx_clav(n) = nanmean(abs(TrialData(n).Results.vCLAV(:,i)));
                elseif i == 2
                    stats.meanVy(n) = nanmean(abs(v));
                    stats.meanPosVy(n) = nanmean(v(indPos));
                    stats.meanNegVy(n) = nanmean(v(indNeg));
                    stats.meanVy_clav(n) = nanmean(abs(TrialData(n).Results.vCLAV(:,i)));
                else
                    stats.meanVz(n) = nanmean(abs(v));
                    stats.meanPosVz(n) = nanmean(v(indPos));
                    stats.meanNegVz(n) = nanmean(v(indNeg));
                    stats.meanVz_clav(n) = nanmean(abs(TrialData(n).Results.vCLAV(:,i)));
                end
                clear indPos indNeg v
            end
            
            %% Metrics for regression of rFin acc, vel, pos to force at interaction point
%             if TrialData(n).Results.px < 0.05 % if fit is sig.
%                 stats.mx(n) = TrialData(n).Results.cx(1);
%                 stats.bx(n) = TrialData(n).Results.cx(2);
%                 stats.kx(n) = TrialData(n).Results.cx(3);
%                 stats.Rsqx(n) = TrialData(n).Results.rsqx;
%             end
%             if TrialData(n).Results.py < 0.05 % if fit is sig.
%                 stats.my(n) = TrialData(n).Results.cy(1);
%                 stats.by(n) = TrialData(n).Results.cy(2);
%                 stats.ky(n) = TrialData(n).Results.cy(3);
%                 stats.Rsqy(n) = TrialData(n).Results.rsqy;
%             end
%             if TrialData(n).Results.pz < 0.05 % if fit is sig.
%                 stats.mz(n) = TrialData(n).Results.cz(1);
%                 stats.bz(n) = TrialData(n).Results.cz(2);
%                 stats.kz(n) = TrialData(n).Results.cz(3);
%                 stats.Rsqz(n) = TrialData(n).Results.rsqz;
%             end
            
            %% Metrics for regression of POB clav acc, vel, pos to force
            if TrialData(n).Results.px_clav < 0.05 % if fit is sig.
                stats.mx_clav(n) = TrialData(n).Results.cx_clav(1);
                stats.bx_clav(n) = TrialData(n).Results.cx_clav(2);
                stats.kx_clav(n) = TrialData(n).Results.cx_clav(3);
                stats.Rsqx_clav(n) = TrialData(n).Results.rsqx_clav;
            end
            if TrialData(n).Results.py_clav < 0.05 % if fit is sig.
                stats.my_clav(n) = TrialData(n).Results.cy_clav(1);
                stats.by_clav(n) = TrialData(n).Results.cy_clav(2);
                stats.ky_clav(n) = TrialData(n).Results.cy_clav(3);
                stats.Rsqy_clav(n) = TrialData(n).Results.rsqy_clav;
            end
            if TrialData(n).Results.pz_clav < 0.05 % if fit is sig.
                stats.mz_clav(n) = TrialData(n).Results.cz_clav(1);
                stats.bz_clav(n) = TrialData(n).Results.cz_clav(2);
                stats.kz_clav(n) = TrialData(n).Results.cz_clav(3);
                stats.Rsqz_clav(n) = TrialData(n).Results.rsqz_clav;
            end
            
            %% Correlation of power to regression residual
            
            % X/ML dir
           stats.rPowerFresX(n) = TrialData(n).Results.corrPowerFresX;
                
            %% Additional metrics that depend on force data
            
            %% Cross-correlation of force to displacement in ML and Vert. dir's
            stats.xcorrX(n) = TrialData(n).Results.xcorrX;
            stats.lagX(n) = TrialData(n).Results.lagX;
            stats.xcorrZ(n) = TrialData(n).Results.xcorrZ;
            stats.lagZ(n) = TrialData(n).Results.lagZ;
            
            % Correlation
            stats.rFvertClav(n) = TrialData(n).Results.rFvertClav; 
            stats.rFlatCOM(n) = TrialData(n).Results.rFlatCOM; 
            stats.rFvertCOM(n) = TrialData(n).Results.rFvertCOM;
            stats.rFlatClav(n) = TrialData(n).Results.rFlatClav;
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
