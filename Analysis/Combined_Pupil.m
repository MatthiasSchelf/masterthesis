%% Load in eeglab
eeglab;

%% Once you have eeglab, run this code. 

% List of participant IDs
participants = {'sub-032'};

% Initialize structures to store participant-specific results
participant_results = struct();

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_'); 
    
    % Load in the data
    pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/%s.mat', participant);
    pupil_data = load(pupil_filename);
    beh_pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/pupil/%s_task-memory_events.tsv', participant, participant);
    beh_pupil = readtable(beh_pupil_filename, 'Delimiter', '\t', 'FileType', 'text');

    % Combine time and diameter into a single matrix
    time = pupil_data.variables_to_export.world_index;
    diameter = pupil_data.variables_to_export.diameter;
    pupil_combined = [time, diameter];
    % Sampling rate
    sampling_rate = 120;
    % Baseline duration
    baseline_duration = 2;

    % Extract Pupil
    [pupillometry_data_blocks_of_five, baseline_blocks_of_five] = extract_pupil_data_baseline_and_interest(beh_pupil, '5', 5, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_nine, baseline_blocks_of_nine] = extract_pupil_data_baseline_and_interest(beh_pupil, '9', 9, pupil_combined, baseline_duration, sampling_rate);
    [pupillometry_data_blocks_of_thirteen, baseline_blocks_of_thirteen] = extract_pupil_data_baseline_and_interest(beh_pupil, '3', 13, pupil_combined, baseline_duration, sampling_rate);

    % Combine Periods of Interest
    average_time_series_pupil = average_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);
    multivariate_time_series_pupil = multivariate_pupil_interest(pupillometry_data_blocks_of_five, pupillometry_data_blocks_of_nine, pupillometry_data_blocks_of_thirteen);

    % Combine Baseline 
    average_baseline_value_pupil = average_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);
    multivariate_baseline_value_pupil = multivariate_pupil_baseline(baseline_blocks_of_five, baseline_blocks_of_nine, baseline_blocks_of_thirteen);
    single_baseline_value_pupil = multivariate_reduce_one(multivariate_baseline_value_pupil);
    
    % Difference Calculation
    multivariate_difference_pupil = multivariate_pupil_difference(multivariate_time_series_pupil, single_baseline_value_pupil);
    difference_pupil = average_time_series_pupil - average_baseline_value_pupil;
    
    % MI, JMI and II average
    MI_pupil_average = calculateMI_pupil(difference_pupil);
    JMI_pupil_average = calculateJMI_pupil(difference_pupil);
    II_pupil_average = calculateII_combined(difference_pupil, difference_pupil);

    % MI, JMI and II multivariate
    MI_pupil_multivariate = calculateMI_pupil(multivariate_difference_pupil);
    JMI_pupil_multivariate = calculateJMI_pupil(multivariate_difference_pupil);
    II_pupil_multivariate = calculateII_combined(multivariate_difference_pupil, multivariate_difference_pupil);
    
    % Store results in the structure
    participant_results.(participant_field_name).MI_pupil_average = MI_pupil_average;
    participant_results.(participant_field_name).JMI_pupil_average = JMI_pupil_average;
    participant_results.(participant_field_name).II_pupil_average = II_pupil_average;
    participant_results.(participant_field_name).MI_pupil_multivariate = MI_pupil_multivariate;
    participant_results.(participant_field_name).JMI_pupil_multivariate = JMI_pupil_multivariate;
    participant_results.(participant_field_name).II_pupil_multivariate = II_pupil_multivariate;

    % Visualize II Pupil 
    visualize_II(II_pupil_multivariate, participant, 'Multivariate');
    visualize_II(II_pupil_average, participant, 'Average');
end

