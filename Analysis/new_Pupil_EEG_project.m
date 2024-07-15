%% Load in eeglab
eeglab;

%% Once EEGlab is loaded, run this code
mydir=pwd;
% List of participant IDs
participants = {'sub-098'};

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    participant_field_name = strrep(participant, '-', '_');

    % EEG data
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

    % Load in the pupil data
    %pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/%s.mat', participant);
    pupil_filename = [mydir,'/Pupilprepro/',participant,'_selected_and_preprocessed.mat'];
    pupil_data = load(pupil_filename);
    %beh_pupil_filename = sprintf('//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/%s/pupil/%s_task-memory_events.tsv', participant, participant);
    beh_pupil_filename = [mydir,'/RawData/',participant,'/pupil/',participant,'_task-memory_events.tsv'];
    beh_pupil = readtable(beh_pupil_filename, 'Delimiter', '\t', 'FileType', 'text');

    % Combine time and pupil dilation
    pupil_time = pupil_data.variables_to_export.pupil_timestamp;
    diameter = pupil_data.variables_to_export.diameter;
    pupil_combined = [pupil_time, diameter];

    % Sampling rate
    sampling_rate = 120;

    % Real code: 
    % Extract relevant time periods for EEG
    [eeg_data_selected_5, eeg_data_selected_9, eeg_data_selected_13] = extract_eeg(beh_eeg, eeg_combined);

    % Extract relevant time periods for Pupil
    [pupil_data_selected_5] = extract_pupil(beh_pupil, '5', 5, pupil_combined, sampling_rate);
    [pupil_data_selected_9] = extract_pupil(beh_pupil, '9', 9, pupil_combined, sampling_rate);
    [pupil_data_selected_13] = extract_pupil(beh_pupil, '3', 13, pupil_combined, sampling_rate);

    % Get the data in the correct format (Trials x Time).
    eeg_data_selected_5_formatted = rearrange_data(eeg_data_selected_5);
    eeg_data_selected_9_formatted = rearrange_data(eeg_data_selected_9);
    eeg_data_selected_13_formatted = rearrange_data(eeg_data_selected_13);
    pupil_data_selected_5_formatted = rearrange_data(pupil_data_selected_5);
    pupil_data_selected_9_formatted = rearrange_data(pupil_data_selected_9);
    pupil_data_selected_13_formatted = rearrange_data(pupil_data_selected_13);

    % Interpolate Nan
    eeg_data_selected_5_formatted = interpolate_nans(eeg_data_selected_5_formatted);
    eeg_data_selected_9_formatted = interpolate_nans(eeg_data_selected_9_formatted);
    eeg_data_selected_13_formatted = interpolate_nans(eeg_data_selected_13_formatted);
    pupil_data_selected_5_formatted = interpolate_nans(pupil_data_selected_5_formatted);
    pupil_data_selected_9_formatted = interpolate_nans(pupil_data_selected_9_formatted);
    pupil_data_selected_13_formatted = interpolate_nans(pupil_data_selected_13_formatted);

    % Trim down to equal lenght
    [eeg_data_selected_5_formatted, pupil_data_selected_5_formatted, eeg_data_selected_9_formatted, pupil_data_selected_9_formatted, eeg_data_selected_13_formatted, pupil_data_selected_13_formatted] = truncate_data(eeg_data_selected_5_formatted, pupil_data_selected_5_formatted, eeg_data_selected_9_formatted, pupil_data_selected_9_formatted, eeg_data_selected_13_formatted, pupil_data_selected_13_formatted);

    % Calculate MI and correlation time course for condition 5 (Trials x
    % Time)
    csddat_5 = eeg_data_selected_5_formatted;
    eyestim_5 = pupil_data_selected_5_formatted;
    [Ntrl, Nt] = size(csddat_5);
    %copula transform for GCMI
    ceeg_5 = copnorm(eeg_data_selected_5_formatted);
    ceyestim_5 = copnorm(pupil_data_selected_5_formatted);

    % repeat calculation at each time point
    rho_5 = zeros(1,Nt);
    I_5 = zeros(1,Nt);
    for ti=1:Nt
        %correlation between pupil  and EEG data at each time point
       rho_5(ti) = corr(csddat_5(:,ti), eyestim_5(:,ti),'type', 'spearman');
        %GCMI between pupil and EEG data at each time point
       I_5(ti) = mi_gg(ceeg_5(:,ti), ceyestim_5(:,ti), true, true);
    end

    % Save results for condition 5
    save_directory = 'D:\UGent_gerelateerd\Masterproef\Data\Output_PP';
    save_filename_rho_5 = [participant '_rho_5.mat'];
    save_filename_I_5 = [participant '_I_5.mat'];
    save_path_rho_5 = fullfile(save_directory, save_filename_rho_5);
    save_path_I_5 = fullfile(save_directory, save_filename_I_5);
    save(save_path_rho_5, 'rho_5');
    save(save_path_I_5, 'I_5');

    % Calculate MI and correlation time course for condition 9 (Trials x
    % Time)
    csddat_9 = eeg_data_selected_9_formatted;
    eyestim_9 = pupil_data_selected_9_formatted;
    [Ntrl, Nt] = size(csddat_9);
    %copula transform for GCMI
    ceeg_9 = copnorm(eeg_data_selected_9_formatted);
    ceyestim_9 = copnorm(pupil_data_selected_9_formatted);

    % repeat calculation at each time point
    rho_9 = zeros(1,Nt);
    I_9 = zeros(1,Nt);
    for ti=1:Nt
        %correlation between pupil stimulus and EEG data at each time point
       rho_9(ti) = corr(csddat_9(:,ti), eyestim_9(:,ti),'type', 'spearman');
        % between pupil s and EEG data at each time point
       I_9(ti) = mi_gg(ceeg_9(:,ti), ceyestim_9(:,ti), true, true);
    end

    % Save results for condition 9
    save_filename_rho_9 = [participant '_rho_9.mat'];
    save_filename_I_9 = [participant '_I_9.mat'];
    save_path_rho_9 = fullfile(save_directory, save_filename_rho_9);
    save_path_I_9 = fullfile(save_directory, save_filename_I_9);
    save(save_path_rho_9, 'rho_9');
    save(save_path_I_9, 'I_9');

    % Calculate MI and correlation time course for condition 13 (Trials x
    % Time)
    csddat_13 = eeg_data_selected_13_formatted;
    eyestim_13 = pupil_data_selected_13_formatted;
    [Ntrl, Nt] = size(csddat_13);
    %copula transform for GCMI
    ceeg_13 = copnorm(eeg_data_selected_13_formatted);
    ceyestim_13 = copnorm(pupil_data_selected_13_formatted);

    % repeat calculation at each time point
    rho_13 = zeros(1,Nt);
    I_13 = zeros(1,Nt);
    for ti=1:Nt
        %correlation between pupil and EEG data at each time point
       rho_13(ti) = corr(csddat_13(:,ti), eyestim_13(:,ti),'type', 'spearman');
        %GCMI between pupil  and EEG data at each time point
       I_13(ti) = mi_gg(ceeg_13(:,ti), ceyestim_13(:,ti), true, true);
    end
    
    % Save results for condition 13
    save_filename_rho_13 = [participant '_rho_13.mat'];
    save_filename_I_13 = [participant '_I_13.mat'];
    save_path_rho_13 = fullfile(save_directory, save_filename_rho_13);
    save_path_I_13 = fullfile(save_directory, save_filename_I_13);
    save(save_path_rho_13, 'rho_13');
    save(save_path_I_13, 'I_13');
end
