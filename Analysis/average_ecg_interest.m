%Functionfor average period of interest of ECG
function average_time_series_ecg = average_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13)
    % Concatenate all ECG data
    all_ecg_data = {ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13};

    % Find the maximum length among all ECG data within each condition
    max_lengths_within_conditions = cellfun(@(x) max(cellfun(@(y) size(y, 1), x)), all_ecg_data);

    % Find the overall maximum length among all conditions
    max_length = max(max_lengths_within_conditions);

    % Initialize a cell array to store padded ECG data for each condition
    padded_data_within_conditions = cell(size(all_ecg_data));

    % Pad the ECG data within each condition to ensure equal length (if
    % unequal = error)
    for i = 1:numel(all_ecg_data)
        data_within_condition = all_ecg_data{i};
        padded_data_within_condition = cell(size(data_within_condition));
        for j = 1:numel(data_within_condition)
            data = data_within_condition{j};
            padded_data_within_condition{j} = [data; NaN(max_length - size(data, 1), 2)];
        end
        padded_data_within_conditions{i} = padded_data_within_condition;
    end

    % Average the padded data across all datasets
    average_time_series_ecg = nanmean(cell2mat(cellfun(@(x) cell2mat(cellfun(@(y) y(:, 2), x, 'UniformOutput', false)), padded_data_within_conditions, 'UniformOutput', false)), 2);
end
