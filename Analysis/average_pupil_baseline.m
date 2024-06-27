function baseline_value_pupil = average_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen)
    % Nested function to calculate the averaged segment
    function average_segment = calculate_averaged_segment(baseline_blocks)
        num_blocks = numel(baseline_blocks);
        if num_blocks == 0
            average_segment = [];
            return;
        end
        
        % Filter out blocks shorter than 15
        baseline_blocks_filtered = baseline_blocks(cellfun(@(x) size(x, 1) >= 15, baseline_blocks));
        num_blocks_filtered = numel(baseline_blocks_filtered);
        
        % If no valid blocks left, return empty
        if num_blocks_filtered == 0
            average_segment = [];
            return;
        end
        
        % Determine the minimum length of remaining segments
        min_length = min(cellfun(@(x) size(x, 1), baseline_blocks_filtered));
        
        % Initialize the averaged segment
        average_segment = zeros(min_length, 1);
        
        % Calculate the averaged time segment
        for i = 1:min_length
            timepoint_values = zeros(num_blocks_filtered, 1);
            for j = 1:num_blocks_filtered
                timepoint_values(j) = baseline_blocks_filtered{j}(i, 2);
            end
            average_segment(i) = mean(timepoint_values);
        end
    end

    % Calculate the averaged time segments for each dataset
    average_baseline_five = calculate_averaged_segment(baseline_blocks_of_five);
    average_baseline_nine = calculate_averaged_segment(baseline_blocks_of_nine);
    average_baseline_thirteen = calculate_averaged_segment(baseline_blocks_of_thirteen);

    % Combine the average segments into a matrix
    average_segments_matrix = [average_baseline_five, ...
                               average_baseline_nine, ...
                               average_baseline_thirteen];

    % Calculate the average segment of the average segments
    average_baseline = mean(average_segments_matrix, 2);

    % Calculate singular average value 
    baseline_value_pupil = mean(average_baseline);
end