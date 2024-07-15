function [eeg_data_5, pupil_data_5, eeg_data_9, pupil_data_9, eeg_data_13, pupil_data_13] = truncate_data(eeg_data_5, pupil_data_5, eeg_data_9, pupil_data_9, eeg_data_13, pupil_data_13)
    % Helper function to truncate data to the minimum size of two matrices
    function [truncated_eeg, truncated_pupil] = truncate_to_min_size(eeg_data, pupil_data)
        min_size = min(size(eeg_data, 2), size(pupil_data, 2));
        truncated_eeg = eeg_data(:, 1:min_size);
        truncated_pupil = pupil_data(:, 1:min_size);
    end

    % Truncate data for each condition
    [eeg_data_5, pupil_data_5] = truncate_to_min_size(eeg_data_5, pupil_data_5);
    [eeg_data_9, pupil_data_9] = truncate_to_min_size(eeg_data_9, pupil_data_9);
    [eeg_data_13, pupil_data_13] = truncate_to_min_size(eeg_data_13, pupil_data_13);
end
