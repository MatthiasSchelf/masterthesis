% Participants
participants = {'sub-098'};
for i = 1:length(participants)
    % Indicate the path where the data is stored.
    data_path = 'D:\UGent_gerelateerd\Masterproef\Data\ECGprepro';

    % Construct the full file path with the correct filename
    participant_id = participants{i};
    file_name = [participant_id, '_task-memory_ecg.mat'];
    file_path = fullfile(data_path, file_name);

    % Load the ECG data
    loaded_data = load(file_path);
    
    % Assume that the loaded data contains a variable 'ecg_signal' and its sampling rate 'fs'
    ecg_signal = loaded_data.HEP.ecg;
    fs = loaded_data.HEP.srate;
    
    % New sampling rate
    new_fs = 120;
    
    % Resample to 120 Hz
    [p, q] = rat(new_fs / fs); % Determine resampling factors
    resampled_ecg_signal = resample(ecg_signal, p, q);
    
    % Save the resampled dataset
    % Create a new file name to indicate that it is resampled
    new_file_name = [participant_id, '_task-memory_ecg_resampled_120.mat'];
    new_file_path = fullfile(data_path, new_file_name);
    
    % Save the resampled ECG signal and the new sampling rate
    save(new_file_path, 'resampled_ecg_signal', 'new_fs');
end
