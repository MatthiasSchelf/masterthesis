%% First, load in EEGlab 
eeglab 

%% Now load in the data 
ECG_data = load("//client/d$/UGent_gerelateerd/Masterproef/Data/ECGprepro/sub-032_task-memory_ecg.mat");

ECG_beh = readtable("//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/sub-032/ecg/sub-032_task-memory_events.tsv", 'Delimiter', '\t', 'FileType', 'text');

%% Extract ECG data
ecg_data = ECG_data.HEP.ecg;

% Create combined matrix
ecg_combined = [(1:length(ecg_data))' ecg_data];

%% Get correct time frames
% Define sequence numbers (5, 9, and 13)
sequence_numbers = [5, 9, 13];

% Initialize cell arrays to store EEG data for each condition
ecg_data_selected_5 = {};
ecg_data_selected_9 = {};
ecg_data_selected_13 = {};

% Loop over sequence numbers
for j = 1:numel(sequence_numbers)
    % Get the current sequence number
    sequence_number = sequence_numbers(j);

    % Identify the onset times and trial types for the current sequence number
    onsets = ECG_beh.onset;
    trial_types = ECG_beh.trial_type;

    % Initialize a cell array to store extracted EEG data for the current condition
    ecg_data_selected = {};

    % Determine time periods and extract EEG data
    is_in_block = false;
    block_start_time = 0;

    for i = 1:numel(onsets)
        onset_time = onsets(i);
        trial_type = trial_types{i};

        % Construct the trial type string based on the current sequence number
        trial_type_string = sprintf('memory 01/%02d', sequence_number);

        % Check if the trial type indicates the start of a block
        if contains(trial_type, trial_type_string) && ~is_in_block
            is_in_block = true;
            block_start_time = onset_time;
        elseif contains(trial_type, sprintf('memory %02d/%02d', sequence_number, sequence_number)) && is_in_block
            % This indicates the end of the block

            % Define the time window for the block (from block start time to onset time)
            window_start = block_start_time;
            window_end = onset_time;

            % Find the indices corresponding to the time window in the EEG data
            indices = ecg_combined(:,1) >= window_start & ecg_combined(:,1) <= window_end;

            % Extract EEG data for this time period
            ecg_data_selected{end+1} = ecg_combined(indices,:);

            % Reset block flag
            is_in_block = false;
        end
    end

    % If a block started but didn't end (i.e., the last block), include all data from the block start time to the end
    if is_in_block
        % Define the time window for the last block (from block start time to end of EEG data)
        window_start = block_start_time;
        window_end = ecg_combined(end, 1);

        % Find the indices corresponding to the time window in the EEG data
        indices = ecg_combined(:,1) >= window_start & ecg_combined(:,1) <= window_end;

        % Extract EEG data for this time period
        ecg_data_selected{end+1} = ecg_combined(indices,:);
    end

    % Store EEG data for the current condition in the appropriate cell array
    switch sequence_number
        case 5
            ecg_data_selected_5 = ecg_data_selected;
        case 9
            ecg_data_selected_9 = ecg_data_selected;
        case 13
            ecg_data_selected_13 = ecg_data_selected;
    end
end

%% Average

% Concatenate all ECG data
all_ecg_data = {ecg_data_selected_5, ecg_data_selected_9, ecg_data_selected_13};

% Find the maximum length among all ECG data within each condition
max_lengths_within_conditions = cellfun(@(x) max(cellfun(@(y) size(y, 1), x(:, 2))), all_ecg_data);

% Find the overall maximum length among all conditions
max_length = max(max_lengths_within_conditions);

% Initialize a cell array to store padded ECG data for each condition
padded_data_within_conditions = cell(size(all_ecg_data));

% Pad the ECG data within each condition to ensure equal length
for i = 1:numel(all_ecg_data)
    data_within_condition = all_ecg_data{i};
    padded_data_within_condition = cell(size(data_within_condition));
    for j = 1:numel(data_within_condition)
        data = data_within_condition{j};
        padded_data_within_condition{j} = [data; NaN(max_length - size(data, 1), 2)];
    end
    padded_data_within_conditions{i} = padded_data_within_condition;
end

% Average the padded data across all datasets
average_time_series_ecg = nanmean(cell2mat(cellfun(@(x) cell2mat(x(:, 2)), padded_data_within_conditions, 'UniformOutput', false)), 2);


%% MI calculations

MI = calculateMI(average_time_series_ecg);

%% JMI calculations

JMI = calculateJMI(average_time_series_ecg);

%% II calculations

II = calculateII(average_time_series_ecg,average_time_series_ecg);

%% Functions 

% Copula normalization function
function cdata = copnorm(data)
    % Copula normalization: converting data to uniform marginals
    [n, m] = size(data);
    cdata = tiedrank(data) / (n + 1); % Empirical CDF ranking
end

% Mutual information calculation function using mi_gg
function MI = calculateMI(data)
    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = size(cdata, 2);

    % Initialize MI array
    MI = zeros(1, ntime_info);

    % Calculate MI for each time point
    for t = 1:ntime_info
        % Extract the time series at time t
        series = cdata(:, t);

        % Check for NaNs or Infs
        if any(isnan(series)) || any(isinf(series))
            warning('Data contains NaNs or Infs at time point %d. Skipping.', t);
            MI(t) = NaN;
            continue;
        end

        % Calculate MI using mi_gg function
        try
            MI(t) = mi_gg(series, series, false); % Placeholder
        catch ME
            warning('MI calculation failed at time point %d: %s', t, ME.message);
            MI(t) = NaN;
        end
    end
end


% Joint mutual information calculation function
function JMI = calculateJMI(data)
    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = size(cdata, 2);

    % Initialize JMI array
    JMI = zeros(ntime_info, ntime_info);

    % Calculate JMI for each pair of time points
    for t1 = 1:ntime_info
        for t2 = (t1 + 1):ntime_info
            % Extract the time series at time t1 and t2
            series1 = cdata(:, t1);
            series2 = cdata(:, t2);

            % Check for NaNs or Infs
            if any(isnan(series1)) || any(isinf(series1)) || any(isnan(series2)) || any(isinf(series2))
                warning('Data contains NaNs or Infs at time points %d and %d. Skipping.', t1, t2);
                JMI(t1, t2) = NaN;
                continue;
            end

            % Calculate JMI using mi_gg function
            try
                JMI(t1, t2) = mi_gg([series1 series2], [series1 series2], false); % Placeholder
            catch ME
                warning('JMI calculation failed at time points %d and %d: %s', t1, t2, ME.message);
                JMI(t1, t2) = NaN;
            end
        end
    end
end

% Define the new function to calculate Interaction Information (II)
function II = calculateII(MI, JMI)
    % Calculate Interaction Information (II) based on MI and JMI
    ntime_info = size(MI, 2);
    II = zeros(ntime_info, ntime_info);
    
    % Loop through pairs of time points
    for t1 = 1:ntime_info
        for t2 = (t1 + 1):ntime_info
            % Calculate II using the formula: II(t1, t2) = JMI(t1, t2) - MI(t1) - MI(t2)
            II(t1, t2) = JMI(t1, t2) - MI(t1) - MI(t2);
        end
    end
    
    % Symmetrize the II matrix
    II = II + II';
end


