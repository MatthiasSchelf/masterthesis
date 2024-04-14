participants = {'sub_098'};
target_frequency = 100;
z_threshold_multiplier = 1;

for i = 1:length(participants)
    % Indicate the path where the data is stored.
    data_path = '\\client\d$\UGent_gerelateerd\Masterproef\Data\Pupilprepro';

    % Construct the full file path with the correct filename
    participant_id = participants{i};
    file_name = [participant_id, '.mat'];
    file_path = fullfile(data_path, file_name);

    % Load the .mat file
    load(file_path);

    % Extract relevant variables from the loaded data
    time = variables_to_export.world_index;
    pupil_size = variables_to_export.diameter;

    % Call the processing function without subsampling
    processed_data = pupil_Preprocessing(time, pupil_size, target_frequency, z_threshold_multiplier);
    
    % Overwrite the original file with the processed data
    variables_to_export.world_index = processed_data.time;
    variables_to_export.diameter = processed_data.pupil_size;
    save(file_path, 'variables_to_export');
end




