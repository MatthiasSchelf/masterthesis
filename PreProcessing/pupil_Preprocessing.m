function pupil_Preprocessing = pupil_Preprocessing(time, pupil_size, target_frequency, z_threshold_multiplier)
    % Input validation
    if isempty(time) || isempty(pupil_size) || numel(time) ~= numel(pupil_size)
        error('Invalid input: time and pupil_size must be non-empty arrays of the same length.');
    end

    % Handle NaN values by interpolating
    non_nan_indices = ~isnan(pupil_size);
    time = time(non_nan_indices);
    pupil_size = pupil_size(non_nan_indices);

    % Downsample the data
    original_frequency = 1 / median(diff(time));
    downsample_factor = round(original_frequency / target_frequency);

    if downsample_factor >= numel(time)
        % If downsampling factor is too large, set it to 1 to keep the original data
        downsample_factor = 1;
    end

    time_to_downsample = time;
    pupil_size_to_downsample = pupil_size;

    time_downsampled = downsample(time_to_downsample, downsample_factor);
    pupil_size_downsampled = downsample(pupil_size_to_downsample, downsample_factor);

    % Dynamic thresholding based on standard deviation
    threshold = z_threshold_multiplier * std(pupil_size_downsampled);

    % Blink detection
    blink_indices = abs(pupil_size_downsampled - mean(pupil_size_downsampled)) > threshold;

    % Median filter for blink removal
    cleaned_pupil_size = medfilt1(pupil_size_downsampled, 5); % Adjust window size as needed

    % Convert time_downsampled and cleaned_pupil_size to double
    time_downsampled = double(time_downsampled);
    cleaned_pupil_size = double(cleaned_pupil_size);

    % Remove duplicate values from time_downsampled and cleaned_pupil_size
    [unique_time, unique_indices] = unique(time_downsampled(~blink_indices));
    unique_pupil_size = cleaned_pupil_size(~blink_indices);
    unique_pupil_size = unique_pupil_size(unique_indices);

    % Interpolation using linear interpolation
    interpolated_pupil_size = interp1(unique_time, unique_pupil_size, time_downsampled, 'linear');

    % Return the function with downsampled time and pupil_size
    pupil_Preprocessing = struct('time', time_downsampled, 'pupil_size', interpolated_pupil_size);
end

