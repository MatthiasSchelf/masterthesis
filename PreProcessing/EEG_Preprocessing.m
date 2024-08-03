function EEG_Preprocessing = EEG_Preprocessing(EEG)
       
    % Step 1: Resample the data
    EEG = pop_resample(EEG, 120);

    % Step 2: High-pass filter at 1 Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'plotfreqz', 1);

    % Step 3: Import channel names
    EEG = pop_chanedit(EEG, 'lookup', 'standard_1005.elc');

    % Step 4: Reject bad channels and indicate the reference channel
    pop_eegplot(EEG, 1, 1, 1);
    EEG = pop_chanedit(EEG, 'append', 66, 'changefield', {67, 'labels', 'FCz'});

    % Step 5: Automatic bad channel rejection using clean_rawdata()
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, 'LineNoiseCriterion', 4, 'Highpass', 'off', 'BurstCriterion', 'off', 'WindowCriterion', 'off');

    % Step 6: Average reference
    EEG = pop_reref(EEG, []);

    % Step 7: ICA
    EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');
    EEG = iclabel(EEG);
    disp(EEG.etc.ic_classification.ICLabel.classes);
    EEG = pop_icflag(EEG, [0 0; 0.8 1; 0.8 1; 0.8 1; 0 0; 0 0; 0 0]);
    
    % Step 8: Get the indices of components marked for rejection
    reject_indices = find(EEG.reject.gcompreject);
    
    % Step 9: Remove rejected components
    if ~isempty(reject_indices)
        EEG.icawinv(:, reject_indices) = 0;  % Set the weights of rejected components to zero
        EEG.icasphere(:, reject_indices) = 0;  % Set the sphere matrix of rejected components to zero
        EEG.icaweights(:, reject_indices) = 0;  % Set the weights of rejected components to zero
        EEG.reject.gcompreject = zeros(1, length(reject_indices));  % Reset the rejection flag
    
        % Save the dataset after component removal
        EEG = pop_saveset(EEG, 'savemode', 'resave');
    else
        disp('No valid components to remove.');
    end

    % Step 10: Interpolate bad channels
    EEG = pop_interp(EEG, find(EEG.reject.rejmanual), 'spherical');

    % Step 11: Save the dataset
    EEG = pop_saveset(EEG, 'savemode', 'resave');

    % Return the processed EEG structure
    EEG_Preprocessing = EEG;
end

