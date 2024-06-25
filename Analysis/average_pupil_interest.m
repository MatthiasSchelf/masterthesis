function average_time_series_pupil = average_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen)

    % Extract pupillometry values from each cell and store in a cell array
    pupil_values_5 = cellfun(@(x) x(:, 2), pupillometry_data_blocks_of_five, 'UniformOutput', false);
    pupil_values_9 = cellfun(@(x) x(:, 2), pupillometry_data_blocks_of_nine, 'UniformOutput', false);
    pupil_values_13 = cellfun(@(x) x(:, 2), pupillometry_data_blocks_of_thirteen, 'UniformOutput', false);

    % Concatenate all pupillometry data
    all_pupil = {pupil_values_5, pupil_values_9, pupil_values_13};

    % Find the length of the longest dataset
    max_length = max(cellfun(@(x) max(cellfun(@(y) length(y), x)), all_pupil));

    % Pad the data with NaN values to ensure equal length
    padded_data = cellfun(@(x) cellfun(@(y) [y; NaN(max_length-length(y), 1)], x, 'UniformOutput', false), all_pupil, 'UniformOutput', false);

    % Concatenate the padded data across all conditions
    concatenated_data = cellfun(@(x) vertcat(x{:}), padded_data, 'UniformOutput', false);

    % Convert concatenated data to a matrix and average across all datasets
    concatenated_matrix = cell2mat(concatenated_data);
    average_time_series = nanmean(concatenated_matrix, 2);

    % Trim the output to the length of the longest column in condition 13
    max_length_13 = max(cellfun(@(x) length(x), pupil_values_13));
    average_time_series_pupil = average_time_series(1:max_length_13);

end
