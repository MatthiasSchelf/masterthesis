% Here we will combine the CMI of EEG Pupil and ECG

%% Load in the data
% Base directory where data files are located
baseDir = '\\client\d$\UGent_gerelateerd\Masterproef\Data\Output_PP_EEG_Pupil_ECG';

% List of participant IDs
participants = {'sub-034', 'sub-036', 'sub-038', 'sub-039', 'sub-040', 'sub-041', 'sub-042', 'sub-043', 'sub-044', 'sub-045', 'sub-046', 'sub-047', 'sub-048', 'sub-049', 'sub-050', 'sub-051', 'sub-053',  'sub-055',  'sub-057', 'sub-058', 'sub-059', 'sub-060',  'sub-062', 'sub-063', 'sub-064', 'sub-065',  'sub-067', 'sub-068', 'sub-069', 'sub-070', 'sub-071',  'sub-073', 'sub-074', 'sub-075', 'sub-076', 'sub-077', 'sub-078', 'sub-079',  'sub-081', 'sub-082',  'sub-084', 'sub-085', 'sub-086', 'sub-087', 'sub-088', 'sub-089', 'sub-091', 'sub-093' 'sub-096', 'sub-097'};

% List of suffixes to loop over
suffixes = {'_5', '_9', '_13'};

% Initialize a structure to store data for each participant
data = struct();

% Loop over each participant
for i = 1:length(participants)
    participant = participants{i};
    
    % Replace invalid characters in participant ID for structure field name
    participantFieldName = strrep(participant, '-', '_');
    
    % Initialize a substructure for the participant
    data.(participantFieldName) = struct();
    
    % Loop over each suffix
    for j = 1:length(suffixes)
        suffix = suffixes{j};
        
        % Construct the full path to the data file
        filePath = fullfile(baseDir, sprintf('%s_CMI%s.mat', participant, suffix));
        
        % Load the data
        loadedData = load(filePath);
        
        % Extract the field names dynamically
        fieldNames = fieldnames(loadedData);
        
        % Store the loaded data directly under the participant's structure
        % without the extra nesting
        data.(participantFieldName).(['CMI' suffix]) = loadedData.(fieldNames{1});
    end
end


% Make sure they all have equal length per condition 

% Ensure that data are equal in length per condition by trimming to the shortest length
for j = 1:length(suffixes)
    suffix = suffixes{j};
    minLength = Inf;
    
    % Find the minimum length for the current suffix
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        % Field name for the current suffix
        fieldName = ['CMI' suffix];
        
        % Update minLength with the length of the current data
        minLength = min([minLength, length(data.(participantFieldName).(fieldName))]);
    end
    
    % Trim data to the minimum length for the current suffix
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        % Field name for the current suffix
        fieldName = ['CMI' suffix];
        
        % Trim the data
        data.(participantFieldName).(fieldName) = data.(participantFieldName).(fieldName)(1:minLength);
    end
end

% Combine all rho and mi per condition 

% Initialize structures to store combined data and averaged results
combined_CMI_5 = struct('mean', [], 'covariance', []);
combined_CMI_9 = struct('mean', [], 'covariance', []);
combined_CMI_13 = struct('mean', [], 'covariance', []);

% Loop over each condition suffix
for j = 1:length(suffixes)
    suffix = suffixes{j};
    
    % Initialize arrays to accumulate data across participants
    all_CMI_data = [];
    
    % Loop over each participant
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        % Get the data for current condition
        CMI_data = data.(participantFieldName).(['CMI' suffix]);
        
        % Concatenate data across participants (transpose to match your data structure)
        all_CMI_data = [all_CMI_data; CMI_data]; % Concatenate row-wise
    end
    
    % Compute the mean vector and covariance matrix across trials
    mean_CMI_data = mean(all_CMI_data, 1); % Mean along rows (participants), result is 1 x time points vector
    cov_CMI_data = cov(all_CMI_data); % Covariance matrix
    
    % Store the mean and covariance in the appropriate variables
    switch suffix
        case '_5'
            combined_CMI_5.mean = mean_CMI_data;
            combined_CMI_5.covariance = cov_CMI_data;
        case '_9'
            combined_CMI_9.mean = mean_CMI_data;
            combined_CMI_9.covariance = cov_CMI_data;
        case '_13'
            combined_CMI_13.mean = mean_CMI_data;
            combined_CMI_13.covariance = cov_CMI_data;
    end
end

% Visualise
% Plot for condition _5
time_points_5 = 1:length(combined_CMI_5.mean);
figure;
plot_mean(time_points_5, combined_CMI_5.mean, 'b');
title('CMI for 5 number sequence');
xlabel('Time Points');
ylabel('CMI');
grid on;

% Plot for condition _9
time_points_9 = 1:length(combined_CMI_9.mean);
figure;
plot_mean(time_points_9, combined_CMI_9.mean, 'g');
title('CMI for 9 number sequence');
xlabel('Time Points');
ylabel('CMI');
grid on;

% Plot for condition _13
time_points_13 = 1:length(combined_CMI_13.mean);
figure;
plot_mean(time_points_13, combined_CMI_13.mean, 'r');
title('CMI for 13 number sequence');
xlabel('Time Points');
ylabel('CMI');
grid on;

% Combined Plot for all conditions
figure;
hold on;
plot_mean(time_points_13, combined_CMI_13.mean, 'r'); % Plot 13 first
plot_mean(time_points_9, combined_CMI_9.mean, 'g');  % Plot 9 second
plot_mean(time_points_5, combined_CMI_5.mean, 'b');  % Plot 5 last
title('CMI for all number sequences');
xlabel('Time Points');
ylabel('CMI');
legend('13 number sequence', '9 number sequence', '5 number sequence');
grid on;
hold off;

% Function to plot mean data only
function plot_mean(x, mean_data, color)
    plot(x, mean_data, color, 'LineWidth', 0.5);
end

