%% First, load in EEGlab 
eeglab;

%% Once EEGlab is loaded, run this code

% List of participant IDs
participants = {'sub-032'};

% Initialize structures to store results
participant_results = struct();

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_'); 
    
    % Load in the data 
    eeg_filename = sprintf('//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/%s_task-memory_eeg.set', participant);
    eeg_data = pop_loadset('filename', eeg_filename);
    beh_eeg_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/eeg/%s_task-memory_events.tsv', participant, participant);
    beh_eeg = readtable(beh_eeg_filename, 'Delimiter', '\t', 'FileType', 'text');
    ecg_filename = sprintf('//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/%s_task-memory_ecg.mat', participant);
    ECG_data = load(ecg_filename);
    beh_ecg_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/ecg/%s_task-memory_events.tsv', participant, participant);
    ECG_beh = readtable(beh_ecg_filename, 'Delimiter', '\t', 'FileType', 'text');
    
    % EEG and ECG preparations
    % Extract EEG data from a specific channel (Fz) along with time
    channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
    if any(channel_index)
        eeg_data_fz = eeg_data.data(channel_index, :);
        num_points = eeg_data.pnts; % Extract number of time points
        sampling_rate = eeg_data.srate; % Extract sampling rate
        eeg_times = (0:num_points-1) / sampling_rate; % Generate time points
        eeg_combined = [eeg_times', eeg_data_fz']; % Combine time points and electrical values
    else
        error('Fz channel not found in the EEG data.');
    end
    % Create matrix with time and ecg
    ecg_data = ECG_data.HEP.ecg;
    ecg_combined = [(1:length(ecg_data))' ecg_data];
    % Define sequence numbers (5, 9, and 13)
    sequence_numbers = [5, 9, 13];
    
    % Extract EEG
    [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13, eeg_baseline_5, eeg_baseline_9, eeg_baseline_13] = extract_eeg_data_baseline_and_interest(beh_eeg, eeg_combined);
    
    % Extract ECG
    [ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13, ecg_baseline_5, ecg_baseline_9, ecg_baseline_13] = extract_ecg_data_baseline_and_interest(sequence_numbers, ECG_beh, ecg_combined);
    
    % Combine periods of interest EEG
    average_time_series_eeg = average_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);
    multivariate_time_series_eeg = multivariate_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);

    % Combine baseline EEG
    average_baseline_value_eeg = average_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);
    multivariate_baseline_value_eeg = multivariate_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);

    % Combine period of interest ECG
    average_time_series_ecg = average_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);
    multivariate_time_series_ecg = multivariate_ecg_interest(ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13);

    % Combine baseline ECG
    average_baseline_value_ecg = average_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);
    multivariate_baseline_value_ecg = multivariate_ecg_baseline(ecg_baseline_5, ecg_baseline_9, ecg_baseline_13);

    % Difference Baseline - Period of Interest
    difference_time_series_eeg = average_time_series_eeg - average_baseline_value_eeg;
    difference_time_series_ecg = average_time_series_ecg - average_baseline_value_ecg;
    multivariate_difference_eeg = multivariate_eeg_difference(multivariate_time_series_eeg, single_baseline_value_eeg);
    multivariate_difference_ecg = multivariate_ecg_difference(multivariate_time_series_ecg, single_baseline_value_ecg);

    % MI JMI and II average
    MI_eeg_average = calculateMI_eeg(difference_time_series_eeg);
    MI_ecg_average = calculateMI_ecg(difference_time_series_ecg);
    JMI_eeg_average = calculateJMI_eeg(difference_time_series_eeg);
    JMI_ecg_average = calculateJMI_ecg(difference_time_series_ecg);
    II_average = calculateII_different(difference_time_series_ecg, difference_time_series_eeg);

    % MI JMI and II multivariate
    MI_ecg_multivariate = calculateMI_ecg(multivariate_difference_ecg);
    JMI_ecg_multivariate = calculateJMI_ecg(multivariate_difference_ecg);
    MI_eeg_multivariate = calculateMI_eeg(multivariate_difference_eeg);
    JMI_eeg_multivariate = calculateJMI_eeg(multivariate_difference_eeg);
    II_multivariate = calculateII_different(multivariate_difference_ecg, multivariate_difference_eeg);

    % Store results in the structure
    participant_results.(participant_field_name).II_eeg_ecg_average = II_average;
    participant_results.(participant_field_name).II_eeg_ecg_multivariate = II_multivariate;

    % Visualize II 
    visualize_II(II_multivariate, participant, 'Multivariate');
    visualize_II(II_average, participant, 'Average');

end
