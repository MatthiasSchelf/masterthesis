%% Load in EEGlab
eeglab;

%% When EEGlab is loaded, run this.

% List of participant IDs
participants = {'sub-032'};

% Initialize structures to store Interactive Information results
participant_results = struct();

% Loop for each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_'); 

    % Load in the data
    pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/%s.mat', participant);
    pupil_data = load(pupil_filename);
    beh_pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/pupil/%s_task-memory_events.tsv', participant, participant);
    beh_pupil = readtable(beh_pupil_filename, 'Delimiter', '\t', 'FileType', 'text');
    eeg_filename = sprintf('//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/%s_task-memory_eeg.set', participant);
    eeg_data = pop_loadset('filename', eeg_filename);
    beh_eeg_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/eeg/%s_task-memory_events.tsv', participant, participant);
    beh_eeg = readtable(beh_eeg_filename, 'Delimiter', '\t', 'FileType', 'text');

    % First take some necessary steps for EEG and Pupil
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
    % Combine Time and Pupillometry
    time = pupil_data.variables_to_export.world_index;
    diameter = pupil_data.variables_to_export.diameter;
    pupil_combined = [time, diameter];
    % Sampling rate
    sampling_rate = 120;
    baseline_duration = 2;

    % Extract EEG
    [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13, eeg_baseline_5, eeg_baseline_9, eeg_baseline_13] = extract_eeg_data_baseline_and_interest(beh_eeg, eeg_combined);

    % Extract Pupil
    [pupillometry_data_blocks_of_five, baseline_blocks_of_five] = extract_pupil_data_baseline_and_interest(beh_pupil, '5', 5, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_nine, baseline_blocks_of_nine] = extract_pupil_data_baseline_and_interest(beh_pupil, '9', 9, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_thirteen, baseline_blocks_of_thirteen] = extract_pupil_data_baseline_and_interest(beh_pupil, '3', 13, pupil_combined, baseline_duration, sampling_rate);

    % Combine Periods of Interest EEG
    average_time_series_eeg = average_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);
    multivariate_time_series_eeg = multivariate_eeg_interest(eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13);

    % Combine Baseline EEG 
    average_baseline_value_eeg = average_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);
    multivariate_baseline_value_eeg = multivariate_eeg_baseline(eeg_baseline_5, eeg_baseline_9, eeg_baseline_13);

    % Combine Periods of Interest Pupil
    average_time_series_pupil = average_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);
    multivariate_time_series_pupil = multivariate_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);

    % Combine Baseline Pupil
    average_baseline_value_pupil = average_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);
    multivariate_baseline_value_pupil = multivariate_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);
  
    % Difference measures 
    difference_time_series_eeg = average_time_series_eeg - average_baseline_value_eeg;
    difference_time_series_pupil = average_time_series_pupil - baseline_value_pupil;
    multivariate_difference_eeg = multivariate_eeg_difference(multivariate_time_series_eeg, single_baseline_value_eeg);
    multivariate_difference_pupil = multivariate_pupil_difference(multivariate_time_series_pupil, single_baseline_value_pupil);

    % MI JMI and II average
    MI_eeg_average = calculateMI_eeg(difference_time_series_eeg);
    MI_pupil_average = calculateMI_pupil(difference_time_series_pupil);
    JMI_eeg_average = calculateJMI_eeg(difference_time_series_eeg);
    JMI_pupil_average = calculateJMI_pupil(difference_time_series_pupil);
    II_average = calculateII_different(difference_time_series_eeg, difference_time_series_pupil);

    % MI JMI and II multivariate
    MI_eeg_multivariate = calculateMI_eeg(multivariate_difference_eeg);
    JMI_eeg_multivariate = calculateJMI_eeg(multivariate_difference_eeg);
    MI_pupil_multivariate = calculateMI_pupil(multivariate_difference_pupil);
    JMI_pupil_multivariate = calculateJMI_pupil(multivariate_difference_pupil);
    II_multivariate = calculateII_different(multivariate_difference_pupil, multivariate_difference_eeg);

    % Store results in the structure
    participant_results.(participant_field_name).II_eeg_pupil_average = II_average;
    participant_results.(participant_field_name).II_eeg_pupil_multivariate = II_multivariate;

    % Visualize II Pupil 
    visualize_II(II_multivariate, participant, 'Multivariate');
    visualize_II(II_average, participant, 'Average');
end
