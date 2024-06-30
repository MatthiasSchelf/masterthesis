%% First, load in EEGlab 
eeglab 

%% Once EEGlab is loaded, run this code.

% List of participant IDs
participants = {'sub-032'};

% Initialize structures to store results
participant_results = struct();

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_'); 

    % Load in the data
    ecg_filename = sprintf('//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/%s_task-memory_ecg.mat', participant);
    ECG_data = load(ecg_filename);
    beh_ecg_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/ecg/%s_task-memory_events.tsv', participant, participant);
    ECG_beh = readtable(beh_ecg_filename, 'Delimiter', '\t', 'FileType', 'text');

    % First necessary steps for ECG
    % Select data column
    ecg_data = ECG_data.HEP.ecg;
    % Create combined matrix
    ecg_combined = [(1:length(ecg_data))' ecg_data];
    % Define sequence numbers (5, 9, and 13)
    sequence_numbers = [5, 9, 13];

    % Extract ECG
    [ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13, ecg_baseline_5, ecg_baseline_9, ecg_baseline_13] = extract_ecg_data_baseline_and_interest(sequence_numbers, ECG_beh, ecg_combined);

    % Combine Periods of Interest ECG
    average_time_series_ecg = average_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);
    multivariate_time_series_ecg = multivariate_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);

    % Combine Baseline ECG
    average_baseline_value_ecg = average_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);
    multivariate_baseline_value_ecg = multivariate_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);
    single_baseline_value_ecg = multivariate_reduce_one(multivariate_baseline_value_ecg);

    % Difference Period of Interest - Baseline
    difference_time_series_ecg = average_time_series_ecg - average_baseline_value_ecg;
    multivariate_difference_ecg = multivariate_ecg_difference(multivariate_time_series_ecg, single_baseline_value_ecg);

    % MI, JMI and II average
    MI_ecg_average = calculateMI_ecg(difference_time_series_ecg);
    JMI_ecg_average = calculateJMI_ecg(difference_time_series_ecg);
    II_ecg_average = calculateII_combined(difference_time_series_ecg, difference_time_series_ecg);
    
    % MI, JMI and II multivariate
    MI_ecg_multivariate = calculateMI_ecg(multivariate_difference_ecg);
    JMI_ecg_multivariate = calculateJMI_ecg(multivariate_difference_ecg);
    II_ecg_multivariate = calculateII_combined(multivariate_difference_ecg, multivariate_difference_ecg);
    
    % Store results in the structure
    participant_results.(participant_field_name).MI_ecg_average = MI_ecg_average;
    participant_results.(participant_field_name).JMI_ecg_average = JMI_ecg_average;
    participant_results.(participant_field_name).II_ecg_average = II_ecg_average;
    participant_results.(participant_field_name).MI_ecg_multivariate = MI_ecg_multivariate;
    participant_results.(participant_field_name).JMI_ecg_multivariate = JMI_ecg_multivariate;
    participant_results.(participant_field_name).II_ecg_multivariate = II_ecg_multivariate;

    % Visualize II ECG
    visualize_II(II_ecg_multivariate, participant, 'Multivariate');
    visualize_II(II_ecg_average, participant, 'Average');
end

