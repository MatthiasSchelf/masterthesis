function [pupillometry_data] = extract_pupil(beh_pupil, digit, blockSize, pupil_combined, sampling_rate)
    % Select and group trials
    beh_pupil.label = string(beh_pupil.label);
    one_before_last_digit = extractBetween(beh_pupil.label, strlength(beh_pupil.label)-1, strlength(beh_pupil.label)-1);
    interesting_rows = beh_pupil(startsWith(beh_pupil.label, '6') & strcmp(one_before_last_digit, digit), :);
    num_interesting_rows = height(interesting_rows);
    blocks = cell(1, ceil(num_interesting_rows / blockSize));
    for i = 1:numel(blocks)
        start_index = (i - 1) * blockSize + 1;
        end_index = min(i * blockSize, num_interesting_rows);
        blocks{i} = interesting_rows(start_index:end_index, :);
    end
    
    % Fill up the trials
    filled_blocks = cellfun(@(block) fillMissingTimestamps(block), blocks, 'UniformOutput', false);
    
    % Extract pupillometry data
    pupillometry_data = cellfun(@(block) extractPupillometryData(block, pupil_combined), filled_blocks, 'UniformOutput', false);

    % Nested function to fill up missing timestamps in a block
    function filled_block = fillMissingTimestamps(block)
        timestamps = block.timestamp;
        start_timestamp = min(timestamps);
        end_timestamp = max(timestamps);
        step_size = 1 / 120;  % More accurate step size calculation
        complete_timestamps = (start_timestamp:step_size:end_timestamp)';
        % Find the missing timestamps
        missing_timestamps = setdiff(complete_timestamps, timestamps);
        % Create rows for the missing timestamps
        num_missing = numel(missing_timestamps);
        missing_rows = array2table(nan(num_missing, width(block)), 'VariableNames', block.Properties.VariableNames);
        missing_rows.timestamp = missing_timestamps;
        % Combine and sort the block
        filled_block = [block; missing_rows];
        filled_block = sortrows(filled_block, 'timestamp');
    end

    % Nested function to extract pupillometry data for each block
    function pupillometry_data = extractPupillometryData(block, pupil_combined)
        timestamps = block.timestamp;
        pupillometry_data = nan(size(block, 1), 2); % Column 1: time, Column 2: diameter
        for i = 1:size(block, 1)
            timestamp = timestamps(i);
            % Use a tolerance for matching timestamps
            block_index = find(abs(pupil_combined(:, 1) - timestamp) < (1 / sampling_rate) / 2, 1); 
            if ~isempty(block_index)
                pupillometry_data(i, 1) = pupil_combined(block_index, 1); % time
                pupillometry_data(i, 2) = pupil_combined(block_index, 2); % diameter
            end
        end
    end
end


