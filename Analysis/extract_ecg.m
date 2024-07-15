function [ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13] = extract_ecg(beh_ecg, ecg_combined)
    % Define sequence numbers (5, 9, and 13)
    sequence_numbers = [5, 9, 13];

    % Initialize cell arrays to store ECG data for each condition
    ecg_data_selected_5 = {};
    ecg_data_selected_9 = {};
    ecg_data_selected_13 = {};

    % Loop over sequence numbers
    for j = 1:numel(sequence_numbers)
        % Get the current sequence number
        sequence_number = sequence_numbers(j);

        % Identify the onset times and trial types for the current sequence number
        onsets = beh_ecg.onset;
        trial_types = beh_ecg.trial_type; % Use beh_ecg.trial_type assuming it contains the labels

        % Initialize cell array to store results
        ecg_data_selected = {};

        % Determine time periods and extract ECG data
        is_in_block = false;
        block_start_time = 0;

        for i = 1:numel(onsets)
            onset_time = onsets(i);
            trial_type = trial_types{i};

            % Construct the trial type string based on the current sequence number
            trial_type_string = sprintf('memory 01/%02d', sequence_number);

            % Check if the trial type indicates the start of a block
            if contains(trial_type, trial_type_string) && ~is_in_block
                is_in_block = true;
                block_start_time = onset_time;
            elseif contains(trial_type, sprintf('memory %02d/%02d', sequence_number, sequence_number)) && is_in_block
                % This indicates the end of the block

                % Define the time window for the block (from block start time to onset time)
                window_start = block_start_time;
                window_end = onset_time;

                % Find the indices corresponding to the time window in the ECG data
                indices = ecg_combined(:, 1) >= window_start & ecg_combined(:, 1) <= window_end;

                % Extract ECG data for this time period
                ecg_data_selected{end+1} = ecg_combined(indices, :);

                % Reset block flag
                is_in_block = false;
            end
        end

        if is_in_block % For the last block 
            % Define the time window for the last block (from block start time to end of ECG data)
            window_start = block_start_time;
            window_end = ecg_combined(end, 1);

            % Find the indices corresponding to the time window in the ECG data
            indices = ecg_combined(:, 1) >= window_start & ecg_combined(:, 1) <= window_end;

            % Extract ECG data for this time period
            ecg_data_selected{end+1} = ecg_combined(indices, :);
        end

        % Store ECG period of interest in the appropriate cell array
        switch sequence_number
            case 5
                ecg_data_selected_5 = ecg_data_selected;
            case 9
                ecg_data_selected_9 = ecg_data_selected;
            case 13
                ecg_data_selected_13 = ecg_data_selected;
        end
    end
end

