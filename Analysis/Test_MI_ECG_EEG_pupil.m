%Load in eeglab 

eeglab 

%% 

% Load in all the necessary datasets.

% Load pupil data from a different directory
load("//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub_032.mat");

% Load ECG data from a different directory
load("//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/sub-032_task-memory_ecg.mat");

% Load EEG data
eeg_data = pop_loadset('filename', '//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/sub-032_task-memory_eeg.set');

%%

% Do the necessary processing 

% For the pupil 

% Extract pupil diameter data
pupil_data = variables_to_export;
time = pupil_data.world_index;
diameter = pupil_data.diameter;

% Combine time and diameter into a single matrix
pupil_combined = [time, diameter];

% For the ECG 

% Extract ECG data from the HEP structure
ecg_data = HEP.ecg;

% For the EEG 

% Extract EEG data from a specific channel (Fz)
channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
if any(channel_index)
    eeg_data_fz = eeg_data.data(channel_index, :);
else
    error('Fz channel not found in the EEG data.');
end

% Check the number of trials (samples) in all three datasets
num_trials_pupil = size(pupil_combined, 1);
num_trials_ecg = size(ecg_data, 1);
num_trials_eeg = size(eeg_data_fz, 2);

% Ensure the number of trials match across all datasets
min_trials = min([num_trials_pupil, num_trials_ecg, num_trials_eeg]);
pupil_combined = pupil_combined(1:min_trials, :);
ecg_data = ecg_data(1:min_trials, :);
eeg_data_fz = eeg_data_fz(:, 1:min_trials);

%Now calculate the MI 

gccmi_ccc(pupil_combined, ecg_data,eeg_data_fz)
