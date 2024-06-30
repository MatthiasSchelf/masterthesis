% This function uses the multivariate probability distribution across
% trials for the periods of interest of ECG
function time_resolved_distributions_ecg = multivariate_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13)
    % Concatenate all ECG data
    all_ecg_data = {ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13};

    % Find the maximum length among all ECG data within each condition
    max_lengths_within_conditions = cellfun(@(x) max(cellfun(@(y) size(y, 1), x)), all_ecg_data);

    % Find the overall maximum length among all conditions
    max_length = max(max_lengths_within_conditions);

    % Initialize a cell array to store padded ECG data for each condition 
    padded_data_within_conditions = cell(size(all_ecg_data));

    % Pad the ECG data within each condition to ensure equal length (no
    % equal length = error).
    for i = 1:numel(all_ecg_data)
        data_within_condition = all_ecg_data{i};
        padded_data_within_condition = cell(size(data_within_condition));
        for j = 1:numel(data_within_condition)
            data = data_within_condition{j};
            padded_data_within_condition{j} = [data; NaN(max_length - size(data, 1), size(data, 2))];
        end
        padded_data_within_conditions{i} = padded_data_within_condition;
    end

    % Concatenate the padded data across all conditions
    concatenated_data_5 = vertcat(padded_data_within_conditions{1}{:});
    concatenated_data_9 = vertcat(padded_data_within_conditions{2}{:});
    concatenated_data_13 = vertcat(padded_data_within_conditions{3}{:});

    % Combine all conditions into one dataset for each time point
    combined_data = cat(3, concatenated_data_5(:, 2), concatenated_data_9(:, 2), concatenated_data_13(:, 2));

    % Get the number of time points and the number of trials
    [num_time_points, ~, num_conditions] = size(combined_data);

    % Initialize a cell array to store results
    time_resolved_distributions_ecg = cell(num_time_points, 1);

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

        % Store the mean and covariance in the previously made cell array
        time_resolved_distributions_ecg{t} = struct('mu', mu, 'Sigma', Sigma);
    end

    % Trim the output to the length of the longest column in condition 13
    % (normally this is automatically OK).
    max_length_13 = max(cellfun(@(x) size(x, 1), ecg_data_selected_13));
    time_resolved_distributions_ecg = time_resolved_distributions_ecg(1:max_length_13);
end
