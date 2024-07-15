% Load EEGlab

eeglab;
%%

% Participants
participants = {'sub-035', 'sub-036',  'sub-038', 'sub-039', 'sub-040', 'sub-041', 'sub-042', 'sub-043', 'sub-044', 'sub-045', 'sub-046', 'sub-047', 'sub-048', 'sub-049', 'sub-050', 'sub-051', 'sub-052', 'sub-053',  'sub-055', 'sub-056', 'sub-057', 'sub-058', 'sub-059', 'sub-060', 'sub-061', 'sub-062', 'sub-063', 'sub-064', 'sub-065', 'sub-067', 'sub-068', 'sub-069', 'sub-070', 'sub-071', 'sub-072', 'sub-073', 'sub-074', 'sub-075', 'sub-076', 'sub-077', 'sub-078', 'sub-079', 'sub-080', 'sub-081', 'sub-082', 'sub-083', 'sub-084', 'sub-085', 'sub-086', 'sub-087', 'sub-088', 'sub-089',  'sub-091', 'sub-092', 'sub-093', 'sub-095', 'sub-096', 'sub-097', 'sub-098'};
for i = 1:length(participants)
    % Indicate the path where the data is stored.
    data_path = 'D:\UGent_gerelateerd\Masterproef\Data\EEGprepro';

    % Construct the full file path with the correct filename
    participant_id = participants{i};
    file_name = [participant_id, '_task-memory_eeg.set'];
    file_path = fullfile(data_path, file_name);

    % Load the data
    EEG = pop_loadset('filename', file_name, 'filepath', data_path);
    
    % Resample to 120 Hz
    EEG = pop_resample(EEG, 120);
    
    % Name the resampled data 
    new_filename = [participants{i} '_task-memory_eeg_processed_120.set'];
    
    % Save the resampled dataset
    pop_saveset(EEG, 'filename', new_filename, 'filepath', data_path);
end
