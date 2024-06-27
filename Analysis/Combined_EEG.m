%% Load in eeglab
eeglab;

%% Once EEGlab is loaded, run this code

% List of participant IDs
participants = {'sub-032'};

% Initialize structures to store participant-specific results
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

    % First take some necessary steps
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

    % Extract EEG
    [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13, eeg_baseline_5, eeg_baseline_9, eeg_baseline_13] = extract_eeg_data_baseline_and_interest(beh_eeg, eeg_combined);

    % Combine Period of interest EEG
    average_time_series_eeg = average_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);
    multivariate_time_series_eeg = multivariate_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);

    % Combine baseline EEG
    average_baseline_value_eeg = average_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);
    multivariate_baseline_value_eeg = multivariate_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);
    single_baseline_value_eeg = multivariate_reduce_one(multivariate_baseline_value_eeg);

    % Calculate Difference
    difference_time_series_eeg = average_time_series_eeg - average_baseline_value_eeg;
    multivariate_difference_eeg = multivariate_eeg_difference(multivariate_time_series_eeg, single_baseline_value_eeg);

    % MI, JMI and II average
    MI_eeg_average = calculateMI_eeg(difference_time_series_eeg);
    JMI_eeg_average = calculateJMI_eeg(difference_time_series_eeg);
    II_eeg_average = calculateII_combined(difference_time_series_eeg, difference_time_series_eeg);

    % MI, JMI and II multivariate
    MI_eeg_multivariate = calculateMI_eeg(multivariate_difference_eeg);
    JMI_eeg_multivariate = calculateJMI_eeg(multivariate_difference_eeg);
    II_eeg_multivariate = calculateII_combined(multivariate_difference_eeg, multivariate_difference_eeg);

    % Store results in the structure
    participant_results.(participant_field_name).MI_eeg_average = MI_eeg_average;
    participant_results.(participant_field_name).JMI_eeg_average = JMI_eeg_average;
    participant_results.(participant_field_name).II_eeg_average = II_eeg_average;
    participant_results.(participant_field_name).MI_eeg_multivariate = MI_eeg_multivariate;
    participant_results.(participant_field_name).JMI_eeg_multivariate = JMI_eeg_multivariate;
    participant_results.(participant_field_name).II_eeg_multivariate = II_eeg_multivariate;

    % Visualize II Pupil 
    visualize_II(II_eeg_multivariate, participant, 'Multivariate');
    visualize_II(II_eeg_average, participant, 'Average');
end
