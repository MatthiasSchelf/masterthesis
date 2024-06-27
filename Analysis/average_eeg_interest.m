function average_time_series_eeg = average_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13)
    % Extract EEG values from each cell and store in a cell array
    eeg_values_5 = cellfun(@(x) x(:, 2), eeg_data_selected_5, 'UniformOutput', false);
    eeg_values_9 = cellfun(@(x) x(:, 2), eeg_data_selected_9, 'UniformOutput', false);
    eeg_values_13 = cellfun(@(x) x(:, 2), eeg_data_selected_13, 'UniformOutput', false);

    % Concatenate all theta power time series data
    all_eeg = {eeg_values_5, eeg_values_9, eeg_values_13};

    % Find the length of the longest dataset
    max_length = max(cellfun(@(x) max(cellfun(@(y) length(y), x)), all_eeg));

    % Pad the data with NaN values to ensure equal length
    padded_data = cellfun(@(x) cellfun(@(y) [y; NaN(max_length-length(y), 1)], x, 'UniformOutput', false), all_eeg, 'UniformOutput', false);

    % Concatenate the padded data across all conditions
    concatenated_data = cellfun(@(x) vertcat(x{:}), padded_data, 'UniformOutput', false);

    % Convert concatenated data to a matrix and average across all datasets
    concatenated_matrix = cell2mat(concatenated_data);
    average_time_series = nanmean(concatenated_matrix, 2);

    % Trim the output to the length of the longest column in condition 13
    max_length_13 = max(cellfun(@(x) length(x), eeg_values_13));
    average_time_series_eeg = average_time_series(1:max_length_13);
end
