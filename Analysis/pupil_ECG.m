%% First load in EEGlab 

eeglab;

%% Once EEGlab is loaded, run this code.

% List of participant IDs
participants = {'sub-032'};

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};

    % Now load in the data
    ecg_filename = sprintf('//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/%s_task-memory_ecg.mat', participant);
    ECG_data = load(ecg_filename);
    beh_ecg_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/ecg/%s_task-memory_events.tsv', participant, participant);
    ECG_beh = readtable(beh_ecg_filename, 'Delimiter', '\t', 'FileType', 'text');
    pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/%s.mat', participant);
    pupil_data = load(pupil_filename);
    beh_pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/pupil/%s_task-memory_events.tsv', participant, participant);
    beh_pupil = readtable(beh_pupil_filename, 'Delimiter', '\t', 'FileType', 'text');

    % Some first steps for ECG and pupil
    % Combine time and ecg
    ecg_data = ECG_data.HEP.ecg;
    ecg_combined = [(1:length(ecg_data))' ecg_data];
    % Combine time and diameter
    time = pupil_data.variables_to_export.world_index;
    diameter = pupil_data.variables_to_export.diameter;
    pupil_combined = [time, diameter];
    % Sampling rate
    sampling_rate = 120;
    % Define sequence numbers (5, 9, and 13)
    sequence_numbers = [5, 9, 13];
    % Baseline duration 
    baseline_duration = 2;

    % Extraction ECG
    [ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13, ecg_baseline_5, ecg_baseline_9, ecg_baseline_13] = extract_ecg_data_baseline_and_interest(sequence_numbers, ECG_beh, ecg_combined);

    % Extraction Pupil
    [pupillometry_data_blocks_of_five, baseline_blocks_of_five] = extract_pupil_data_baseline_and_interest(beh_pupil, '5', 5, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_nine, baseline_blocks_of_nine] = extract_pupil_data_baseline_and_interest(beh_pupil, '9', 9, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_thirteen, baseline_blocks_of_thirteen] = extract_pupil_data_baseline_and_interest(beh_pupil, '3', 13, pupil_combined, baseline_duration, sampling_rate);

    % Combine Period of Interest ECG
    average_time_series_ecg = average_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);
    multivariate_time_series_ecg = multivariate_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);

    % Combine baseline ECG
    average_baseline_value_ecg = average_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);
    multivariate_baseline_value_ecg = multivariate_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);

    % Combine period of interest Pupil
    average_time_series_pupil = average_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);
    multivariate_time_series_pupil = multivariate_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);

    % Combine baseline Pupil
    average_baseline_value_pupil = average_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);
    multivariate_baseline_value_pupil = multivariate_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);

    % Get the difference measures 
    difference_time_series_ecg = average_time_series_ecg - average_baseline_value_ecg;
    difference_time_series_pupil = average_time_series_pupil - baseline_value_pupil;
    multivariate_difference_pupil = multivariate_pupil_difference(multivariate_time_series_pupil, single_baseline_value_pupil);
    multivariate_difference_ecg = multivariate_ecg_difference(multivariate_time_series_ecg, single_baseline_value_ecg);
    
    % MI, JMI, II for average
    MI_ecg_average = calculateMI_ecg(difference_time_series_ecg);
    MI_pupil_average = calculateMI_pupil(difference_time_series_pupil);
    JMI_ecg_average = calculateJMI_ecg(difference_time_series_ecg);
    JMI_pupil_average = calculateJMI_pupil(difference_time_series_pupil);
    II_average = calculateII(difference_time_series_ecg, difference_time_series_pupil);

    % MI, JMI, II for multivariate
    MI_pupil_multivariate = calculateMI_pupil(multivariate_difference_pupil);
    JMI_pupil_multivariate = calculateJMI_pupil(multivariate_difference_pupil);
    MI_ecg_multivariate = calculateMI_ecg(multivariate_difference_ecg);
    JMI_ecg_multivariate = calculateJMI_ecg(multivariate_difference_ecg);
    II_multivariate = calculateII(multivariate_difference_pupil, multivariate_difference_ecg);
end
