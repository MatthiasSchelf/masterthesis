function time_resolved_distributions = multivariate_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13)
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

    % Combine all conditions into one dataset for each time point
    combined_baseline_data = cat(3, concatenated_baseline_5, concatenated_baseline_9, concatenated_baseline_13);

    % Get the number of time points and the number of trials
    [num_time_points, num_trials] = size(combined_baseline_data(:,:,1));

    % Initialize a cell array to hold the distributions for each time point
    time_resolved_distributions = cell(num_time_points, 1);

    % Iterate through each time point to calculate the multivariate distribution
    for t = 1:num_time_points
        % Collect data across trials and conditions for the current time point
        data_at_time_t = [];
        for c = 1:size(combined_baseline_data, 3)
            data_at_time_t = [data_at_time_t; combined_baseline_data(t,:,c)'];
        end

        % Calculate the mean and covariance matrix for the multivariate distribution
        mu = mean(data_at_time_t, 1);
        Sigma = cov(data_at_time_t);

        % Store the mean and covariance in the cell array
        time_resolved_distributions{t} = struct('mu', mu, 'Sigma', Sigma);
    end
end

