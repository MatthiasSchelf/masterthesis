function time_resolved_distributions = multivariate_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13)
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
    concatenated_data_5 = vertcat(padded_data{1}{:});
    concatenated_data_9 = vertcat(padded_data{2}{:});
    concatenated_data_13 = vertcat(padded_data{3}{:});

    % Combine all conditions into one dataset for each time point
    combined_data = cat(3, concatenated_data_5, concatenated_data_9, concatenated_data_13);

    % Get the number of time points and the number of trials
    [num_time_points, ~, num_conditions] = size(combined_data);

    % Initialize a cell array to hold the distributions for each time point
    time_resolved_distributions = cell(num_time_points, 1);

    % Iterate through each time point to calculate the multivariate distribution
    for t = 1:num_time_points
        % Collect data across trials for the current time point
        data_at_time_t = squeeze(combined_data(t, :, :));

        % Remove NaN values
        data_at_time_t = data_at_time_t(~any(isnan(data_at_time_t), 2), :);

        % Calculate the mean and covariance matrix for the multivariate distribution
        if ~isempty(data_at_time_t)
            mu = mean(data_at_time_t, 1);
            Sigma = cov(data_at_time_t);
        else
            mu = NaN(1, num_conditions);
            Sigma = NaN(num_conditions);
        end

        % Store the mean and covariance in the cell array
        time_resolved_distributions{t} = struct('mu', mu, 'Sigma', Sigma);
    end

    % Trim the output to the length of the longest column in condition 13
    max_length_13 = max(cellfun(@(x) length(x), eeg_values_13));
    time_resolved_distributions = time_resolved_distributions(1:max_length_13);
end


