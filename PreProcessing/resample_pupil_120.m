participants = {'sub-032'};

for i = 1:length(participants)
    % Indicate the path where the data is stored.
    data_path = 'D:\UGent_gerelateerd\Masterproef\Data\Pupilprepro';

    % Construct the full file path with the correct filename
    participant_id = participants{i};
    file_name = [participant_id, '_resampled.mat'];
    file_path = fullfile(data_path, file_name);

    % Load the .mat file
    load(file_path);

    % Define the target frequency
    target_frequency = 120; % Hz

    % Get the original sampling rate and time vector
    original_time = processed_data.time;
    original_diameter = processed_data.pupil_size;
    
    % Calculate the original sampling rate
    original_sampling_rate = 1 / mean(diff(original_time));

    % Calculate the new time vector
    new_time = original_time(1):(1/target_frequency):original_time(end);

    % Resample the pupil size data
    resampled_diameter = resample(original_diameter, target_frequency, original_sampling_rate);

    % Create a structure to hold the resampled data
    variables_to_export.world_index = new_time;
    variables_to_export.diameter = resampled_diameter;

    % Save the resampled data back to the .mat file
    save(file_path, 'variables_to_export');
end

