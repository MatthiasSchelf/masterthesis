% Function to extract baseline and period of interest for Pupillometry. 
function [pupillometry_data, baseline_blocks] = extract_pupil_data_baseline_and_interest(beh_pupil, digit, blockSize, pupil_combined, baseline_duration, sampling_rate)
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
    
    % Fill up the trials (sometime only timepoint 4 and 7 are given because
    % their a number is given to the participant, however timepoints 5 and
    % 6 are also of importance since we want to see the complete period and
    % not only the selected time fragments in which a number is given)
    filled_blocks = cellfun(@(block) fillMissingIndices(block), blocks, 'UniformOutput', false);
    
    % Extract pupillometry data
    pupillometry_data = cellfun(@(block) extractPupillometryData(block, pupil_combined), filled_blocks, 'UniformOutput', false);
    
    % Extract baseline data
    baseline_length = baseline_duration * sampling_rate; % Baseline length in timepoints
    baseline_blocks = cell(size(pupillometry_data));
    
    for i = 1:numel(pupillometry_data)
        block = pupillometry_data{i};
        if ~isempty(block)
            period_start_time = block(1, 1);
            start_time = period_start_time - baseline_length;
            if start_time < 1
                start_time = 1;
            end
            end_time = period_start_time - 1;
            baseline_data = pupil_combined(start_time:end_time, :);
            % Ensure unique pairs of timepoints and pupillometry data
            [~, unique_indices] = unique(baseline_data, 'rows', 'stable');
            baseline_blocks{i} = baseline_data(unique_indices, :);
        else
            baseline_blocks{i} = [];
        end
    end
    
    % Nested function to fill up missing indices in a block (again same as
    % previously, all timepoints need to be included).
    function filled_block = fillMissingIndices(block)
        indices = block{:, 1};
        start_index = indices(1);
        end_index = indices(end);
        complete_indices = (start_index:end_index)';
        missing_indices = setdiff(complete_indices, indices);
        missing_rows = array2table([missing_indices, nan(numel(missing_indices), width(block)-1)], 'VariableNames', block.Properties.VariableNames);
        filled_block = [block; missing_rows];
        filled_block = sortrows(filled_block);
    end

    % Nested function to extract pupillometry data for each block
    function pupillometry_data = extractPupillometryData(block, pupil_combined)
        indices = block{:, 1};
        pupillometry_data = zeros(size(block, 1), 2); % Column 1: time, Column 2: diameter
        for i = 1:size(block, 1)
            index = indices(i);
            block_index = find(pupil_combined(:, 1) == index, 1);
            if ~isempty(block_index)
                pupillometry_data(i, 1) = pupil_combined(block_index, 1); % time
                pupillometry_data(i, 2) = pupil_combined(block_index, 2); % diameter
            else
                pupillometry_data(i, :) = NaN; % handle missing data
            end
        end
    end
end