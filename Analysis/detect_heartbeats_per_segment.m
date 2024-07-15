function num_heartbeats_per_segment = detect_heartbeats_per_segment(ecg_data_selected, sampling_rate)
    % Peak detection parameters (adjust as per your signal characteristics)
    threshold = 0.3; % Adjust threshold as per your data
    
    num_segments = numel(ecg_data_selected);
    num_heartbeats_per_segment = zeros(1, num_segments);

    % Loop through each ECG segment and detect heartbeats
    for seg_idx = 1:num_segments
        ecg_segment = ecg_data_selected{seg_idx}(:, 2); % Extract ECG signal from the second column
        [~, peak_locs] = findpeaks(ecg_segment, 'MinPeakHeight', threshold, 'MinPeakDistance', ceil(0.6 * sampling_rate));
        num_heartbeats_per_segment(seg_idx) = length(peak_locs);
    end
end
