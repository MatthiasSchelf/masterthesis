function ecg_matrix = rearrange_ecg(ecg_data_selected)
    % Ensure ecg_data_selected is a cell array
    if ~iscell(ecg_data_selected)
        ecg_data_selected = num2cell(ecg_data_selected);
    end

    num_trials = length(ecg_data_selected);
    % Determine the maximum number of heartbeats (should be 1 for single values)
    max_heartbeats = 1;
    % Initialize the matrix with NaNs
    ecg_matrix = NaN(num_trials, max_heartbeats);
    
    for trial_idx = 1:num_trials
        % Get the ECG data for the current trial
        ecg_trial_data = ecg_data_selected{trial_idx};
        % Assign the ECG data to the corresponding row
        ecg_matrix(trial_idx, 1) = ecg_trial_data;
    end
end

