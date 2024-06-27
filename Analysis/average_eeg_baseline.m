function average_baseline_value_eeg = average_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13)
    % Extract EEG values from each cell and store in a cell array
    baseline_values_5 = cellfun(@(x) x(:, 2), eeg_baseline_5, 'UniformOutput', false);
    baseline_values_9 = cellfun(@(x) x(:, 2), eeg_baseline_9, 'UniformOutput', false);
    baseline_values_13 = cellfun(@(x) x(:, 2), eeg_baseline_13, 'UniformOutput', false);

    % Concatenate all baseline data
    all_baseline = {baseline_values_5, baseline_values_9, baseline_values_13};

    % Since all baseline periods have the same length, we do not need to pad with NaN
    % Concatenate the baseline data across all conditions
    concatenated_baseline_5 = horzcat(baseline_values_5{:});
    concatenated_baseline_9 = horzcat(baseline_values_9{:});
    concatenated_baseline_13 = horzcat(baseline_values_13{:});

    % Calculate the average baseline for each condition
    average_baseline_5 = mean(concatenated_baseline_5, 2);
    average_baseline_9 = mean(concatenated_baseline_9, 2);
    average_baseline_13 = mean(concatenated_baseline_13, 2);

    % Average all baselines together to get a single average baseline period
    average_baseline_series_eeg = mean([average_baseline_5, average_baseline_9, average_baseline_13], 2);

    % Calculate the average of the average baseline series as one number
    average_baseline_value_eeg = mean(average_baseline_series_eeg);
end
