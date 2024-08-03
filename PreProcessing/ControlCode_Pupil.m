participants = {'sub-035', 'sub-036',  'sub-038', 'sub-039', 'sub-040', 'sub-041', 'sub-042', 'sub-043', 'sub-044', 'sub-045', 'sub-046', 'sub-047', 'sub-048', 'sub-049', 'sub-050', 'sub-051', 'sub-052', 'sub-053',  'sub-055', 'sub-056', 'sub-057', 'sub-058', 'sub-059', 'sub-060', 'sub-061', 'sub-062', 'sub-063', 'sub-064', 'sub-065', 'sub-067', 'sub-068', 'sub-069', 'sub-070', 'sub-071', 'sub-072', 'sub-073', 'sub-074', 'sub-075', 'sub-076', 'sub-077', 'sub-078', 'sub-079', 'sub-080', 'sub-081', 'sub-082', 'sub-083', 'sub-084', 'sub-085', 'sub-086', 'sub-087', 'sub-088', 'sub-089',  'sub-091', 'sub-092', 'sub-093', 'sub-095', 'sub-096', 'sub-097', 'sub-098'};
target_frequency = 120;
z_threshold_multiplier = 1;

for i = 1:length(participants)
    % Indicate the path where the data is stored.
    data_path = 'D:\UGent_gerelateerd\Masterproef\Data\Pupilprepro';

    % Construct the full file path with the correct filename
    participant_id = participants{i};
    file_name = [participant_id, '_selected.mat'];
    file_path = fullfile(data_path, file_name);

    % Load the .mat file
    load(file_path);

    % Extract relevant variables from the loaded data
    time = variables_to_export.pupil_timestamp;
    pupil_size = variables_to_export.diameter;

    % Call the processing function without subsampling
    processed_data = pupil_Preprocessing(time, pupil_size, target_frequency, z_threshold_multiplier);
    
    % Prepare the new file name for saving the processed data
    new_file_name = [participant_id, '_selected_and_preprocessed.mat'];
    new_file_path = fullfile(data_path, new_file_name);

    % Save the processed data into a new file
    variables_to_export.pupil_timestamp = processed_data.time;
    variables_to_export.diameter = processed_data.pupil_size;
    save(new_file_path, 'variables_to_export');
end




