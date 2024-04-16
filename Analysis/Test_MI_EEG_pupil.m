%Load in EEGlab

eeglab

%%

% Load EEG data
eeg_data = pop_loadset('filename', '//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/sub-032_task-memory_eeg.set');

% Load pupil data
load("//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub_032.mat");

%%

% Extract pupil diameter data
pupil_data = variables_to_export;
time = pupil_data.world_index;
diameter = pupil_data.diameter;

% Combine time and diameter into a single matrix
pupil_combined = [time, diameter];

% Extract EEG data from a specific channel (Fz)
channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
if any(channel_index)
    eeg_data_fz = eeg_data.data(channel_index, :);
else
    error('Fz channel not found in the EEG data.');
end

% Check the number of trials (samples) in both EEG data and pupil data
num_trials_x = size(eeg_data_fz, 2); % Assuming each column is a trial/sample
num_trials_y = size(pupil_combined, 1);

% If the number of trials doesn't match, adjust the data accordingly
if num_trials_x ~= num_trials_y
    % Adjust the data to match the number of trials
    min_trials = min(num_trials_x, num_trials_y);
    eeg_data_fz = eeg_data_fz(:, 1:min_trials);
    pupil_combined = pupil_combined(1:min_trials, :);
end

% Now you can use the gcmi_cc function
gcmi_score = gcmi_cc(pupil_combined, eeg_data_fz);
disp(['GCMi Score: ', num2str(gcmi_score)]);

