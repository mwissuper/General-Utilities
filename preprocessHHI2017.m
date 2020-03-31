function [TrialData] = preprocessHHI2017(TrialData)
    %PREPROCESSHHI2017: Reorganizes Trial Data from 2017 HHI Experiment
    %   preprocessHHI2017 performs the following actions on each of the
    %   elements in the TrialData array (each of the individual trials):
    %       1) Removes 'Unknown' from Marker Data
    %       2) Reorganizes EMG to contain POB (person on beam) and AP (assistance provider) fields, each with
    %       BCEP, TCEP, WFLE, and WEXT subfields. EMG no longer contains
    %       the ElectricPotential field
    %       3) Pulls Fx, Fy, Fz, Mx, My, Mz from Other.ElectricPotential
    %       and puts them in a new field called Forces, together with
    %       Other.SampleRate. The Other field is removed. 
    %       4) Relabels specific subjects as generic. POB_xx, AP_xx, and
    %       FH_xx are relabled as POB, AP, and FH, respectively. 
    
    %   Luke Drnach
    %   Neuromechanics Laboratory
    %   Georgia Tech and Emory
          
    forcenames = {'Fx','Fy','Fz','Mx','My','Mz'};
    for n = 1:length(TrialData)
        if isfield(TrialData(n),'EMG') && isfield(TrialData(n).EMG,'Voltage')
            % First: Create new fields in the EMG data field
            TrialData(n).EMG.POB = struct();
            TrialData(n).EMG.AP = struct();
            emg_names = fieldnames(TrialData(n).EMG.Voltage);
            for m = 1:length(emg_names)
                if emg_names{m}(end) == 'L'
                    TrialData(n).EMG.POB.(emg_names{m}(1:end-2)) = TrialData(n).EMG.Voltage.(emg_names{m});
                elseif emg_names{m}(end) == 'R'
                    TrialData(n).EMG.AP.(emg_names{m}(1:end-2)) = TrialData(n).EMG.Voltage.(emg_names{m});
                end
            end
            % Second: remove ElectricPotential and Voltage from the EMG field
            TrialData(n).EMG = rmfield(TrialData(n).EMG,{'Voltage','ElectricPotential'});
        end
        % Third: Create a new field for the forces and torques
        if isfield(TrialData(n),'Other') && isfield(TrialData(n).Other,'ElectricPotential')
            TrialData(n).Forces = struct();
            for m = 1:length(forcenames)
                TrialData(n).Forces.(forcenames{m}) = TrialData(n).Other.ElectricPotential.(forcenames{m});
            end
            TrialData(n).Forces.samplerate = TrialData(n).Other.samplerate;
        end
        % Fourth: Check if there are any Unknowns in the Markers field, and
        % remove them
        if isfield(TrialData(n),'Markers')
            MarkerFields = fieldnames(TrialData(n).Markers);
            idx = strncmp(MarkerFields,'Unknown',7);
            if any(idx)
                TrialData(n).Markers = rmfield(TrialData(n).Markers,MarkerFields(idx));
            end           
            
        end
        
    end
    if isfield(TrialData,'Other')
        %Finally, remove the Other field from every trial
        TrialData = rmfield(TrialData,'Other');
    end
    
    % Now, relabel all the subjects with generic names
    ExpNumber = TrialData(1).Info.Subject_1;
    idx = find(ExpNumber == '_');
    ExpNumber = ExpNumber(idx:end);
    oldnames = {['AP',ExpNumber],['POB',ExpNumber],['FH',ExpNumber]};
    newnames = {'AP','POB','FH'};
    field = 'Markers';
    TrialData = relabelData(TrialData,field,oldnames,newnames);
    
end

