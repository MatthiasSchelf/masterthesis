%% Load in eeglab
eeglab;

%% Once EEGlab is loaded, run this code
mydir = pwd;
% List of participant IDs
participants = {'sub-039'};

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_');
    
    % ECG DATA
    ecg_filename = [mydir,'/ECGprepro/',participant,'_task-memory_ecg_resampled_120.mat'];
    ecg_data = load(ecg_filename);
    
    % Extract the resampled ECG signal
    resampled_ecg_signal = ecg_data.resampled_ecg_signal;
    
    % Create time vector starting from 0 with steps of 1/120 (assuming 120 Hz sampling rate)
    num_samples = length(resampled_ecg_signal);
    time_vector = (0:num_samples-1) * (1/120);
    
    % Ensure that time_vector is a column vector and resampled_ecg_signal is also a column vector
    time_vector = time_vector(:); % Ensure it's a column vector
    resampled_ecg_signal = resampled_ecg_signal(:); % Ensure it's a column vector
    
    % Combine time vector and resampled_ecg_signal into ecg_combined
    ecg_combined = [time_vector, resampled_ecg_signal];
    
    beh_ecg_filename = [mydir,'/RawData/',participant,'/ecg/',participant,'_task-memory_events.tsv'];
    beh_ecg = readtable(beh_ecg_filename, 'Delimiter', '\t', 'FileType', 'text');
    
    % EEG DATA
    eeg_filename = [mydir,'/EEGprepro/',participant,'_task-memory_eeg_processed_120.set'];
    eeg_data = pop_loadset('filename', eeg_filename);
    beh_eeg_filename = [mydir,'/RawData/',participant,'/eeg/',participant,'_task-memory_events.tsv'];
    beh_eeg = readtable(beh_eeg_filename, 'Delimiter', '\t', 'FileType', 'text');
    
    % First take some necessary steps
    % Extract EEG data from a specific channel (Fz) along with time
    channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
    if any(channel_index)
        eeg_data_fz = double(eeg_data.data(channel_index, :));
        num_points = eeg_data.pnts; % Extract number of time points
        sampling_rate = 120;
        eeg_times = (0:num_points-1) / sampling_rate; % Generate time points
        eeg_combined = [eeg_times', eeg_data_fz']; % Combine time points and electrical values
    else
        error('Fz channel not found in the EEG data.');
    end
    
    % PUPIL DATA
    pupil_filename = [mydir,'/Pupilprepro/',participant,'_selected_and_preprocessed.mat'];
    pupil_data = load(pupil_filename);
    beh_pupil_filename = [mydir,'/RawData/',participant,'/pupil/',participant,'_task-memory_events.tsv'];
    beh_pupil = readtable(beh_pupil_filename, 'Delimiter', '\t', 'FileType', 'text');
    
    % Combine time and pupil dilation
    pupil_time = pupil_data.variables_to_export.pupil_timestamp;
    diameter = pupil_data.variables_to_export.diameter;
    pupil_combined = [pupil_time, diameter];
    
    % Real code
    % Extract relevant time periods for ECG
    [ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13] = extract_ecg(beh_ecg, ecg_combined);
    
    % Extract relevant time periods for EEG
    [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13] = extract_eeg(beh_eeg, eeg_combined);
    
    % Extract relevant time periods for Pupil
    [pupil_data_selected_5] = extract_pupil(beh_pupil, '5', 5, pupil_combined, sampling_rate);
    [pupil_data_selected_9] = extract_pupil(beh_pupil, '9', 9, pupil_combined, sampling_rate);
    [pupil_data_selected_13] = extract_pupil(beh_pupil, '3', 13, pupil_combined, sampling_rate);
    
    % Get amount of heartbeats per segment
    sampling_rate = 120;
    num_heartbeats_5 = detect_heartbeats_per_segment(ecg_data_selected_5, sampling_rate);
    num_heartbeats_9 = detect_heartbeats_per_segment(ecg_data_selected_9, sampling_rate);
    num_heartbeats_13 = detect_heartbeats_per_segment(ecg_data_selected_13, sampling_rate);
    
    % Get the data in the correct format (Trials x Time).
    eeg_data_selected_5_formatted = rearrange_data(eeg_data_selected_5);
    eeg_data_selected_9_formatted = rearrange_data(eeg_data_selected_9);
    eeg_data_selected_13_formatted = rearrange_data(eeg_data_selected_13);
    pupil_data_selected_5_formatted = rearrange_data(pupil_data_selected_5);
    pupil_data_selected_9_formatted = rearrange_data(pupil_data_selected_9);
    pupil_data_selected_13_formatted = rearrange_data(pupil_data_selected_13);
    ecg_data_selected_5_formatted = rearrange_ecg(num_heartbeats_5);
    ecg_data_selected_9_formatted = rearrange_ecg(num_heartbeats_9);
    ecg_data_selected_13_formatted = rearrange_ecg(num_heartbeats_13);
    
    % Interpolate Nan
    eeg_data_selected_5_formatted = interpolate_nans(eeg_data_selected_5_formatted);
    eeg_data_selected_9_formatted = interpolate_nans(eeg_data_selected_9_formatted);
    eeg_data_selected_13_formatted = interpolate_nans(eeg_data_selected_13_formatted);
    pupil_data_selected_5_formatted = interpolate_nans(pupil_data_selected_5_formatted);
    pupil_data_selected_9_formatted = interpolate_nans(pupil_data_selected_9_formatted);
    pupil_data_selected_13_formatted = interpolate_nans(pupil_data_selected_13_formatted);
    
    % Trim down eeg and pupil to equal length
    [eeg_data_selected_5_formatted, pupil_data_selected_5_formatted, eeg_data_selected_9_formatted, pupil_data_selected_9_formatted, eeg_data_selected_13_formatted, pupil_data_selected_13_formatted] = truncate_data(eeg_data_selected_5_formatted, pupil_data_selected_5_formatted, eeg_data_selected_9_formatted, pupil_data_selected_9_formatted, eeg_data_selected_13_formatted, pupil_data_selected_13_formatted);
    
    % Define save directory
    save_directory = 'D:\UGent_gerelateerd\Masterproef\Data\Output_PP_EEG_PUPIL_ECG'; % replace with actual save directory path
    
    % Conditional mutual information (CMI) for 5
    % copula transform for GCMI
    ceeg_5 = copnorm(eeg_data_selected_5_formatted);
    ceyestim_5 = copnorm(pupil_data_selected_5_formatted);
    [Ntrl, Nt] = size(ceeg_5);
    
    % Now integrate the heartbeat to filter it out.
    CMI_5 = zeros(1,Nt);
    for ti=1:Nt
        % cmi function is like mi but takes three copula normalised inputs
        CMI_5(ti) = cmi_ggg(ceeg_5(:,ti), ceyestim_5(:,ti), ecg_data_selected_5_formatted(:,1), true, true);
    end
    
    % Save results for condition 5
    save_filename_CMI_5 = [participant '_CMI_5.mat'];
    save_path_CMI_5 = fullfile(save_directory, save_filename_CMI_5);
    save(save_path_CMI_5, 'CMI_5');
    
    % Conditional mutual information (CMI) for 9
    % copula transform for GCMI
    ceeg_9 = copnorm(eeg_data_selected_9_formatted);
    ceyestim_9 = copnorm(pupil_data_selected_9_formatted);
    [Ntrl, Nt] = size(ceeg_9);
    
    % Now integrate the heartbeat to filter it out.
    CMI_9 = zeros(1,Nt);
    for ti=1:Nt
        % cmi function is like mi but takes three copula normalised inputs
        CMI_9(ti) = cmi_ggg(ceeg_9(:,ti), ceyestim_9(:,ti), ecg_data_selected_9_formatted(:,1), true, true);
    end
    
    % Save results for condition 9
    save_filename_CMI_9 = [participant '_CMI_9.mat'];
    save_path_CMI_9 = fullfile(save_directory, save_filename_CMI_9);
    save(save_path_CMI_9, 'CMI_9');
    
    % Conditional mutual information (CMI) for 13
    % copula transform for GCMI
    ceeg_13 = copnorm(eeg_data_selected_13_formatted);
    ceyestim_13 = copnorm(pupil_data_selected_13_formatted);
    [Ntrl, Nt] = size(ceeg_13);
    
    % Now integrate the heartbeat to filter it out.
    CMI_13 = zeros(1,Nt);
    for ti=1:Nt
        % cmi function is like mi but takes three copula normalised inputs
        CMI_13(ti) = cmi_ggg(ceeg_13(:,ti), ceyestim_13(:,ti), ecg_data_selected_13_formatted(:,1), true, true);
    end
    
    % Save results for condition 13
    save_filename_CMI_13 = [participant '_CMI_13.mat'];
    save_path_CMI_13 = fullfile(save_directory, save_filename_CMI_13);
    save(save_path_CMI_13, 'CMI_13');
end

