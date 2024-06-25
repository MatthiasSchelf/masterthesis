function time_resolved_distributions_pupil_baseline = multivariate_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen)
    % Nested function to filter and interpolate segments
    function interpolated_segments = filter_and_interpolate_segments(baseline_blocks)
        % Filter out blocks shorter than 15
        baseline_blocks_filtered = baseline_blocks(cellfun(@(x) size(x, 1) >= 15, baseline_blocks));
        
        % If no valid blocks left, return empty
        if isempty(baseline_blocks_filtered)
            interpolated_segments = [];
            return;
        end
        
        % Interpolate segments to a fixed length of 16
        fixed_length = 16;
        interpolated_segments = cellfun(@(x) interp1(1:size(x, 1), x(:, 2), linspace(1, size(x, 1), fixed_length), 'linear', 'extrap'), baseline_blocks_filtered, 'UniformOutput', false);
    end

    % Process each set of baseline blocks
    interpolated_segments_5 = filter_and_interpolate_segments(baseline_blocks_of_five);
    interpolated_segments_9 = filter_and_interpolate_segments(baseline_blocks_of_nine);
    interpolated_segments_13 = filter_and_interpolate_segments(baseline_blocks_of_thirteen);

    % If any condition has no valid segments, initialize empty output
    if isempty(interpolated_segments_5) || isempty(interpolated_segments_9) || isempty(interpolated_segments_13)
        time_resolved_distributions_pupil_baseline = [];
        return;
    end

    % Convert interpolated segments to matrices
    concatenated_data_5 = vertcat(interpolated_segments_5{:});
    concatenated_data_9 = vertcat(interpolated_segments_9{:});
    concatenated_data_13 = vertcat(interpolated_segments_13{:});

    % Ensure that all concatenated data have the same number of trials
    min_trials = min([size(concatenated_data_5, 1), size(concatenated_data_9, 1), size(concatenated_data_13, 1)]);
    concatenated_data_5 = concatenated_data_5(1:min_trials, :);
    concatenated_data_9 = concatenated_data_9(1:min_trials, :);
    concatenated_data_13 = concatenated_data_13(1:min_trials, :);

    % Combine all conditions into one dataset for each time point
    combined_data = cat(3, concatenated_data_5, concatenated_data_9, concatenated_data_13);

    % Get the number of time points and the number of trials
    [num_trials, num_time_points] = size(combined_data(:, :, 1));

    % Initialize a cell array to hold the distributions for each time point
    time_resolved_distributions_pupil_baseline = cell(num_time_points, 1);

    % Iterate through each time point to calculate the multivariate distribution
    for t = 1:num_time_points
        % Collect data across trials and conditions for the current time point
        data_at_time_t = [];
        for c = 1:size(combined_data, 3)
            data_at_time_t = [data_at_time_t; combined_data(:, t, c)];
        end

        % Calculate the mean and covariance matrix for the multivariate distribution
        mu = mean(data_at_time_t, 1);
        Sigma = cov(data_at_time_t);

        % Store the mean and covariance in the cell array
        time_resolved_distributions_pupil_baseline{t} = struct('mu', mu, 'Sigma', Sigma);
    end
end
