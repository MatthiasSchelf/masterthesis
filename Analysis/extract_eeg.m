% Function to extract period of interest for EEG
function [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13] = extract_eeg(beh_eeg, eeg_combined)
    % Define sequence numbers (5, 9, and 13)
    sequence_numbers = [5, 9, 13];

    % Initialize cell arrays to store EEG data for each condition
    eeg_data_selected_5 = {};
    eeg_data_selected_9 = {};
    eeg_data_selected_13 = {};

    % Loop over sequence numbers
    for j = 1:numel(sequence_numbers)
        % Get the current sequence number
        sequence_number = sequence_numbers(j);

        % Identify the onset times and trial types for the current sequence number
        onsets = beh_eeg.onset;
        trial_types = beh_eeg.trial_type;

        % Initialize cell array to store results
        eeg_data_selected = {};

        % Determine time periods and extract EEG data
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

                % Find the indices corresponding to the time window in the EEG data
                indices = eeg_combined(:, 1) >= window_start & eeg_combined(:, 1) <= window_end;

                % Extract EEG data for this time period
                eeg_data_selected{end+1} = eeg_combined(indices, :);

                % Reset block flag
                is_in_block = false;
            end
        end

        if is_in_block % For the last block 
            % Define the time window for the last block (from block start time to end of EEG data)
            window_start = block_start_time;
            window_end = eeg_combined(end, 1);

            % Find the indices corresponding to the time window in the EEG data
            indices = eeg_combined(:, 1) >= window_start & eeg_combined(:, 1) <= window_end;

            % Extract EEG data for this time period
            eeg_data_selected{end+1} = eeg_combined(indices, :);
        end

        % Correctly store EEG period of interest
        switch sequence_number
            case 5
                eeg_data_selected_5 = eeg_data_selected;
            case 9
                eeg_data_selected_9 = eeg_data_selected;
            case 13
                eeg_data_selected_13 = eeg_data_selected;
        end
    end
end
