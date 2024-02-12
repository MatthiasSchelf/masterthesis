%% load eeglab 

eeglab

%% Loop over all the datasets. 

% Common part of the dataset path
common_path = 'D:/UGent_gerelateerd/Masterproef/Data/RawData/';

% Variable parts for the last two blocks of the path
participants = {'sub-032/eeg/sub-032_task-rest_eeg'};

% Loop over participants
for i = 1:length(participants)
   
    % Construct the full dataset path
    data_path = fullfile(common_path, participants{i});

    % Get a list of all .set files in the directory
    set_files = dir(fullfile(data_path, '*.set'));

    % Loop over each .set file
    for j = 1:length(set_files)
        % Load EEG data using EEGLAB functions
        EEG = pop_loadset('filename', set_files(j).name, 'filepath', data_path);

        % Call the processing function
        processed_data_test = EEG_Preprocessing(EEG);

        % Save the processed_data as needed
        save_results(processed_data_test, data_path);
    end
end


