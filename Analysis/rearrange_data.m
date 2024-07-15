function trials_x_time_matrix = rearrange_data(eeg_data_selected)
    num_trials = length(eeg_data_selected);
    % Determine the maximum number of time points
    max_timepoints = max(cellfun(@(x) size(x, 1), eeg_data_selected));
    % Initialize the matrix with NaNs
    trials_x_time_matrix = NaN(num_trials, max_timepoints);
    
    for trial_idx = 1:num_trials
        % Get the EEG data for the current trial
        eeg_trial_data = eeg_data_selected{trial_idx}(:, 2)';
        % Determine the length of the current trial
        trial_length = length(eeg_trial_data);
        % Assign the EEG data to the corresponding row, padded with NaNs if necessary
        trials_x_time_matrix(trial_idx, 1:trial_length) = eeg_trial_data;
    end
end