%% load eeglab + open restEEG

eeglab

%% Loop over all the participants in memory. 

% Variable parts for the last two blocks of the path
participants = {'sub-078','sub-079','sub-080','sub-081','sub-082','sub-083','sub-084','sub-085','sub-086','sub-087','sub-088','sub-089','sub-091','sub-092','sub-093','sub-095','sub-096','sub-097','sub-098'};

% Loop over participants
for i = 1:length(participants)
    
    data_path = 'D:\UGent_gerelateerd\Masterproef\Data\EEGprepro';

    % Construct the full file path
    participant_id = participants{i};
    file_path = fullfile(data_path, [participant_id '_task-memory_eeg.set']);

    % Load EEG data using EEGLAB functions
    EEG = pop_loadset(file_path);

    % Call the processing function
    processed_data = EEG_Preprocessing(EEG);
end



