% Here we will combine the MI of the searate participants into one.

%% Load in the data
% Base directory where data files are located
baseDir = '\\client\d$\UGent_gerelateerd\Masterproef\Data\Output_PP_EEG_Pupil';

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
        
        % Construct the full paths to the data files for 'I' and 'rho'
        filePath_I = fullfile(baseDir, sprintf('%s_I%s.mat', participant, suffix));
        filePath_rho = fullfile(baseDir, sprintf('%s_rho%s.mat', participant, suffix));
        
        % Load the 'I' and 'rho' data
        loadedData_I = load(filePath_I);
        loadedData_rho = load(filePath_rho);
        
        % Extract the field names dynamically
        fieldNames_I = fieldnames(loadedData_I);
        fieldNames_rho = fieldnames(loadedData_rho);
        
        % Store the loaded data directly under the participant's structure
        % without the extra nesting
        data.(participantFieldName).(['I' suffix]) = loadedData_I.(fieldNames_I{1});
        data.(participantFieldName).(['rho' suffix]) = loadedData_rho.(fieldNames_rho{1});
    end
end

% Make sure they all have equal length per condition 

% Ensure that I and rho data are equal in length per condition by trimming to the shortest length
for j = 1:length(suffixes)
    suffix = suffixes{j};
    minLength = Inf;
    
    % Find the minimum length for I and rho data for the current suffix
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        I_field = ['I' suffix];
        rho_field = ['rho' suffix];
        
        minLength = min([minLength, length(data.(participantFieldName).(I_field)), length(data.(participantFieldName).(rho_field))]);
    end
    
    % Trim I and rho data to the minimum length for the current suffix
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        I_field = ['I' suffix];
        rho_field = ['rho' suffix];
        
        % Trim I data
        data.(participantFieldName).(I_field) = data.(participantFieldName).(I_field)(1:minLength);
        
        % Trim rho data
        data.(participantFieldName).(rho_field) = data.(participantFieldName).(rho_field)(1:minLength);
    end
end
% Combine all rho and mi per condition 

% Initialize structures to store combined data and averaged results
combined_I_5 = [];
combined_I_9 = [];
combined_I_13 = [];
combined_rho_5 = [];
combined_rho_9 = [];
combined_rho_13 = [];

% Loop over each condition suffix
for j = 1:length(suffixes)
    suffix = suffixes{j};
    
    % Initialize arrays to accumulate data across participants
    all_I_data = [];
    all_rho_data = [];
    
    % Loop over each participant
    for i = 1:length(participants)
        participant = participants{i};
        participantFieldName = strrep(participant, '-', '_');
        
        % Get the data for current condition and type
        I_data = data.(participantFieldName).(['I' suffix]);
        rho_data = data.(participantFieldName).(['rho' suffix]);
        
        % Concatenate data across participants (transpose to match your data structure)
        all_I_data = [all_I_data; I_data]; % Concatenate row-wise
        all_rho_data = [all_rho_data; rho_data]; % Concatenate row-wise
    end
    
    % Compute the mean vector and covariance matrix across trials
    mean_I_data = mean(all_I_data, 1); % Mean along rows (participants), result is 1 x time points vector
    mean_rho_data = mean(all_rho_data, 1); % Mean along rows (participants), result is 1 x time points vector
    
    % Compute the covariance matrix
    cov_I_data = cov(all_I_data); % Covariance matrix
    cov_rho_data = cov(all_rho_data); % Covariance matrix
    
    % Store the mean and covariance in the appropriate variables
    switch suffix
        case '_5'
            combined_I_5.mean = mean_I_data;
            combined_I_5.covariance = cov_I_data;
            combined_rho_5.mean = mean_rho_data;
            combined_rho_5.covariance = cov_rho_data;
        case '_9'
            combined_I_9.mean = mean_I_data;
            combined_I_9.covariance = cov_I_data;
            combined_rho_9.mean = mean_rho_data;
            combined_rho_9.covariance = cov_rho_data;
        case '_13'
            combined_I_13.mean = mean_I_data;
            combined_I_13.covariance = cov_I_data;
            combined_rho_13.mean = mean_rho_data;
            combined_rho_13.covariance = cov_rho_data;
    end
end

% Visualise


% Define the sampling rate and convert time points to milliseconds
sampling_rate = 120; % in Hz
time_points_5 = 1:length(combined_I_5.mean);
time_in_ms_5 = (time_points_5 - 1) * (1000 / sampling_rate); % Convert to milliseconds

% Plot for condition _5 - I
% Visualise
figure;
plot_mean(time_in_ms_5, combined_I_5.mean, 'b');
title('MI for 5 number sequence');
xlabel('Time (ms)');
ylabel('MI');
grid on;

% Plot for condition _5 - rho
figure;
plot_mean(time_in_ms_5, combined_rho_5.mean, 'b');
title('Spearman Correlation for 5 number sequence');
xlabel('Time (ms)');
ylabel('Correlation');
grid on;


% Define the sampling rate and convert time points to milliseconds
sampling_rate = 120; % in Hz
time_points_9 = 1:length(combined_I_9.mean);
time_in_ms_9 = (time_points_9 - 1) * (1000 / sampling_rate); % Convert to milliseconds

% Plot for condition _9 - I
figure;
plot_mean(time_in_ms_9, combined_I_9.mean, 'g');
title('MI for 9 number sequence');
xlabel('Time (ms)');
ylabel('MI');
grid on;

% Plot for condition _9 - rho
figure;
plot_mean(time_in_ms_9, combined_rho_9.mean, 'g');
title('Spearman Correlation for 9 number sequence');
xlabel('Time (ms)');
ylabel('Correlation');
grid on;


% Define the sampling rate and convert time points to milliseconds
sampling_rate = 120; % in Hz
time_points_13 = 1:length(combined_I_13.mean);
time_in_ms_13 = (time_points_13 - 1) * (1000 / sampling_rate); % Convert to milliseconds

% Plot for condition _13 - I
figure;
plot_mean(time_in_ms_13, combined_I_13.mean, 'r');
title('MI for 13 number sequence');
xlabel('Time (ms)');
ylabel('MI');
grid on;

% Plot for condition _13 - rho
figure;
plot_mean(time_in_ms_13, combined_rho_13.mean, 'r');
title('Spearman Correlation for 13 number sequence');
xlabel('Time (ms)');
ylabel('Correlation');
grid on;

% Visualize
% Plot for MI with all lengths
figure;
hold on;
plot_mean(time_in_ms_13, combined_I_13.mean, 'r'); % Plot 13 first
plot_mean(time_in_ms_9, combined_I_9.mean, 'g');   % Plot 9 second
plot_mean(time_in_ms_5, combined_I_5.mean, 'b');   % Plot 5 last
title('MI Overlap');
xlabel('Time (ms)');
ylabel('MI');
legend({'13 number sequence', '9 number sequence', '5 number sequence'}, 'Location', 'Best');
grid on;
hold off;

% Plot for Spearman Correlation with all lengths
figure;
hold on;
plot_mean(time_in_ms_13, combined_rho_13.mean, 'r'); % Plot 13 first
plot_mean(time_in_ms_9, combined_rho_9.mean, 'g');   % Plot 9 second
plot_mean(time_in_ms_5, combined_rho_5.mean, 'b');   % Plot 5 last
title('Spearman Correlation Overlap');
xlabel('Time (ms)');
ylabel('Correlation');
legend({'13 number sequence', '9 number sequence', '5 number sequence'}, 'Location', 'Best');
grid on;
hold off;

% Function to plot mean data only
function plot_mean(x, mean_data, color)
    plot(x, mean_data, color, 'LineWidth', 0.5);
end

