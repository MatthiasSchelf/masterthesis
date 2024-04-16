
%Load in EEGlab

eeglab

%%

% Load ECG data
load("//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/sub-032_task-memory_ecg.mat");

% Load EEG data
eeg_data = pop_loadset('filename', '//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/sub-032_task-memory_eeg.set');

%%

% Extract EEG data from a specific channel (Fz)
channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
if any(channel_index)
    eeg_data_fz = eeg_data.data(channel_index, :);
else
    error('Fz channel not found in the EEG data.');
end

% Extract ECG data from the HEP structure
ecg_data = HEP.ecg;

% Ensure ECG and EEG data have the same time axis
% Assuming you have a common time axis in both ECG and EEG data
% Adjust timestamps or indices to align the data

% Example: if you have timestamps in milliseconds
ecg_timestamps = 0:1:length(ecg_data)-1; % Assuming ECG data is sampled at 1 Hz
eeg_timestamps = (0:1:length(eeg_data_fz)-1) * (1000 / eeg_data.srate); % Convert EEG indices to milliseconds

% Interpolate ECG data to match EEG timestamps
ecg_data_aligned = interp1(ecg_timestamps, ecg_data, eeg_timestamps, 'linear', 'extrap');

% Now you can use the gcmi_cc function
gcmi_score = gcmi_cc(ecg_data_aligned, eeg_data_fz);
disp(['GCMi Score: ', num2str(gcmi_score)]);
