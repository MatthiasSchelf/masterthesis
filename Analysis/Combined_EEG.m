%% Load in eeglab
eeglab

%% Load in the data

% Load EEG data
eeg_data = pop_loadset('filename', '//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/sub-032_task-memory_eeg.set');

% Load behavior file for eeg
beh_eeg = readtable("//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/sub-032/eeg/sub-032_task-memory_events.tsv", 'Delimiter', '\t', 'FileType', 'text');

%% First take some necessary steps

% Extract EEG data from a specific channel (Fz) along with time
channel_index = strcmpi({eeg_data.chanlocs.labels}, 'Fz');
if any(channel_index)
    eeg_data_fz = eeg_data.data(channel_index, :);
    num_points = eeg_data.pnts; % Extract number of time points
    sampling_rate = eeg_data.srate; % Extract sampling rate
    eeg_times = (0:num_points-1) / sampling_rate; % Generate time points
    eeg_combined = [eeg_times', eeg_data_fz']; % Combine time points and electrical values
else
    error('Fz channel not found in the EEG data.');
end

%% Get correct time frames
% Define sequence numbers (5, 9, and 13)
sequence_numbers = [5, 9, 13];

% Initialize cell arrays to store EEG data for each condition
eeg_data_selected_5 = {};
eeg_data_selected_9 = {};
eeg_data_selected_13 = {};

% Loop over sequence numbers
for j = 1:numel(sequence_numbers)
    % Get the current sequence number
    sequence_number = sequence_numbers(j);

    % Identify the onset times and trial types for the current sequence number
    onsets = beh_eeg.onset;
    trial_types = beh_eeg.trial_type;

    % Initialize a cell array to store extracted EEG data for the current condition
    eeg_data_selected = {};

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
            indices = eeg_combined(:,1) >= window_start & eeg_combined(:,1) <= window_end;

            % Extract EEG data for this time period
            eeg_data_selected{end+1} = eeg_combined(indices,:);

            % Reset block flag
            is_in_block = false;
        end
    end

    % If a block started but didn't end (i.e., the last block), include all data from the block start time to the end
    if is_in_block
        % Define the time window for the last block (from block start time to end of EEG data)
        window_start = block_start_time;
        window_end = eeg_combined(end, 1);

        % Find the indices corresponding to the time window in the EEG data
        indices = eeg_combined(:,1) >= window_start & eeg_combined(:,1) <= window_end;

        % Extract EEG data for this time period
        eeg_data_selected{end+1} = eeg_combined(indices,:);
    end

    % Store EEG data for the current condition in the appropriate cell array
    switch sequence_number
        case 5
            eeg_data_selected_5 = eeg_data_selected;
        case 9
            eeg_data_selected_9 = eeg_data_selected;
        case 13
            eeg_data_selected_13 = eeg_data_selected;
    end
end

%% Get theta in blocks

% Define parameters for time window (e.g., duration and overlap)
window_duration = 0.5; % Duration of each time window in seconds
overlap = 0.25; % Overlap between consecutive time windows as a fraction of window duration
theta_band = [4;8];

% Calculate theta power time series for condition 5
theta_power_time_series_5 = calculate_theta_power_time_series(eeg_data_selected_5, theta_band, sampling_rate, window_duration, overlap);

% Calculate theta power time series for condition 9
theta_power_time_series_9 = calculate_theta_power_time_series(eeg_data_selected_9, theta_band, sampling_rate, window_duration, overlap);

% Calculate theta power time series for condition 13
theta_power_time_series_13 = calculate_theta_power_time_series(eeg_data_selected_13, theta_band, sampling_rate, window_duration, overlap);

%% Average 

% Concatenate all theta power time series data
all_theta_power = {theta_power_time_series_5, theta_power_time_series_9, theta_power_time_series_13};

% Extract EEG values from each cell and store in a cell array
eeg_values = cellfun(@(x) cell2mat(x(:, 2)), all_theta_power, 'UniformOutput', false);

% Find the length of the longest dataset
max_length = max(cellfun(@(x) size(x, 1), eeg_values));

% Pad the data with NaN values to ensure equal length
padded_data = cellfun(@(x) [x; NaN(max_length-size(x, 1), 1)], eeg_values, 'UniformOutput', false);

% Average the padded data across all datasets
average_time_series = nanmean(cell2mat(padded_data), 2);

%% MI calculations

MI = calculateMI(average_time_series);

%% JMI calculations

JMI = calculateJMI(average_time_series);

%% II calculations

II = calculateII(average_time_series,average_time_series);

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



%Function to calculate theta power time series for a given condition
function theta_power_time_series = calculate_theta_power_time_series(eeg_data_selected, theta_band, sampling_rate, window_duration, overlap)
    % Initialize cell array to store theta power time series for each block
    theta_power_time_series = cell(size(eeg_data_selected));

    % Calculate power in theta band for each time window within each block
    for i = 1:numel(eeg_data_selected)
        % Extract EEG data for the block
        eeg_block = eeg_data_selected{i};

        % Extract timepoints and EEG data
        timepoints = eeg_block(:, 1);
        eeg_data = eeg_block(:, 2);

        % Define time windows
        window_start = timepoints(1);
        window_end = window_start + window_duration;

        theta_power_time_series_block = []; % Initialize array to store theta power time series

        % Calculate power in theta band for each time window
        while window_end <= timepoints(end)
            % Find indices corresponding to current time window
            window_indices = timepoints >= window_start & timepoints < window_end;

            % Extract EEG data for current time window
            eeg_window = eeg_data(window_indices);

            % Perform FFT to convert EEG data to frequency domain
            fft_data = fft(eeg_window);

            % Calculate power spectral density
            psd = abs(fft_data).^2;

            % Find frequencies corresponding to FFT bins
            freq = (0:length(psd)-1) * (sampling_rate / length(psd));

            % Find indices corresponding to theta band frequencies
            theta_indices = freq >= theta_band(1) & freq <= theta_band(2);

            % Calculate power in theta band for current time window
            theta_power = sum(psd(theta_indices));

            % Store theta power in time series array
            theta_power_time_series_block = [theta_power_time_series_block; theta_power];

            % Move to next time window
            window_start = window_start + window_duration * (1 - overlap);
            window_end = window_end + window_duration * (1 - overlap);
        end

        % Store theta power time series for current block
        theta_power_time_series{i} = theta_power_time_series_block;
    end
end

