%% load eeglab + open restEEG

eeglab

%% Loop over all the participants in rest. 

% Variable parts for the last two blocks of the path
participants = {'sub-032'};

% Loop over participants
for i = 1:length(participants)
    
    data_path = '\\client\d$\UGent_gerelateerd\Masterproef\Data\EEGprepro';

    % Construct the full file path
    participant_id = participants{i};
    file_path = fullfile(data_path, [participant_id '_task-memory_eeg.set']);

    % Load EEG data using EEGLAB functions
    EEG = pop_loadset(file_path);

    % Call the processing function
    processed_data = EEG_Preprocessing(EEG);
end



