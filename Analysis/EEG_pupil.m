%% Load in eeglab
eeglab

%% Load in the data
pupil_data = load("//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub_032.mat");
beh_pupil = readtable("//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/sub-032/pupil/sub-032_task-memory_events.tsv", 'Delimiter', '\t', 'FileType', 'text');
eeg_data = pop_loadset('filename', '//client/d$/UGent_gerelateerd/Masterproef/Data/EEGprepro/sub-032_task-memory_eeg.set');
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

time = pupil_data.variables_to_export.world_index;
diameter = pupil_data.variables_to_export.diameter;

% Combine time and diameter into a single matrix
pupil_combined = [time, diameter];

% Define sequence numbers (5, 9, and 13)
sequence_numbers = [5, 9, 13];

%% Selecting periods of interest

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

blocks_of_five = selectAndGroupTrials(beh_pupil, '5', 5);
blocks_of_nine = selectAndGroupTrials(beh_pupil, '9', 9);
blocks_of_thirteen = selectAndGroupTrials(beh_pupil, '3', 13);

filled_blocks_of_five = cellfun(@fillMissingIndices, blocks_of_five, 'UniformOutput', false);
filled_blocks_of_nine = cellfun(@fillMissingIndices, blocks_of_nine, 'UniformOutput', false);
filled_blocks_of_thirteen = cellfun(@fillMissingIndices, blocks_of_thirteen, 'UniformOutput', false);

pupillometry_data_blocks_of_five = cellfun(@(block) extractPupillometryData(block, pupil_combined), filled_blocks_of_five, 'UniformOutput', false);
pupillometry_data_blocks_of_nine = cellfun(@(block) extractPupillometryData(block, pupil_combined), filled_blocks_of_nine, 'UniformOutput', false);
pupillometry_data_blocks_of_thirteen = cellfun(@(block) extractPupillometryData(block, pupil_combined), filled_blocks_of_thirteen, 'UniformOutput', false);

%% Averaging

% Concatenate all pupillometry data
all_pupillometry_data = [pupillometry_data_blocks_of_five(:); pupillometry_data_blocks_of_nine(:); pupillometry_data_blocks_of_thirteen(:)];

% Extract pupillometry values from each cell and store in a cell array
pupil_values = cellfun(@(x) cell2mat(x), all_pupillometry_data, 'UniformOutput', false);

% Find the length of the longest dataset
max_length = max(cellfun(@(x) size(x, 1), pupil_values));

% Pad the data with NaN values to ensure equal length
padded_data = cellfun(@(x) [x; NaN(max_length-size(x, 1), size(x, 2))], pupil_values, 'UniformOutput', false);

% Average the padded data across all datasets
average_time_series_pupil = nanmean(cat(3, padded_data{:}), 3);

% Concatenate all theta power time series data
all_theta_power = {theta_power_time_series_5, theta_power_time_series_9, theta_power_time_series_13};

% Extract EEG values from each cell and store in a cell array
eeg_values = cellfun(@(x) cell2mat(x(:, 2)), all_theta_power, 'UniformOutput', false);

% Find the length of the longest dataset
max_length = max(cellfun(@(x) size(x, 1), eeg_values));

% Pad the data with NaN values to ensure equal length
padded_data = cellfun(@(x) [x; NaN(max_length-size(x, 1), 1)], eeg_values, 'UniformOutput', false);

% Average the padded data across all datasets
average_time_series_eeg = nanmean(cell2mat(padded_data), 2);

%% MI  JMI and II

MI_eeg = calculateMI(average_time_series_eeg);
JMI_eeg = calculateJMI(average_time_series_eeg);

MI_pupil = calculateMI(average_time_series_pupil);
JMI_pupil = calculateJMI(average_time_series_pupil);

II = calculateII(average_time_series_eeg,average_time_series_pupil);

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



% Define a function to select and group trials
function blocks = selectAndGroupTrials(beh_pupil, digit, blockSize)
    % Convert 'label' column to text
    beh_pupil.label = string(beh_pupil.label);

    % Extracting the one-before-last digit from the "label" column
    one_before_last_digit = extractBetween(beh_pupil.label, strlength(beh_pupil.label)-1, strlength(beh_pupil.label)-1);

    % Extracting rows where the third column "label" starts with '6' and the one-before-last digit matches the specified digit
    interesting_rows = beh_pupil(startsWith(beh_pupil.label, '6') & strcmp(one_before_last_digit, digit), :);

    % Group the interesting rows in blocks
    num_interesting_rows = height(interesting_rows);
    blocks = cell(1, ceil(num_interesting_rows / blockSize));
    for i = 1:numel(blocks)
        start_index = (i - 1) * blockSize + 1;
        end_index = min(i * blockSize, num_interesting_rows);
        blocks{i} = interesting_rows(start_index:end_index, :);
    end
end

% Define a function to fill up missing indices in a block
function filled_block = fillMissingIndices(block)
    % Get the indices from the first column
    indices = block{:, 1};

    % Determine the range of indices
    start_index = indices(1);
    end_index = indices(end);

    % Generate a complete set of indices
    complete_indices = (start_index:end_index)';

    % Find the missing indices
    missing_indices = setdiff(complete_indices, indices);

    % Create rows with missing indices and empty values for other columns
    missing_rows = array2table([missing_indices, nan(numel(missing_indices), width(block)-1)], 'VariableNames', block.Properties.VariableNames);

    % Combine the original block with the rows with missing indices
    filled_block = [block; missing_rows];

    % Sort the filled block based on the indices
    filled_block = sortrows(filled_block);
end

% Extract pupillometry data for each block and condition
function pupillometry_data = extractPupillometryData(block, pupil_combined)
    % Get the indices from the first column
    indices = block{:, 1};
    
    % Initialize pupillometry data
    pupillometry_data = cell(size(block, 1), 1);
    
    % Loop through each row in the block
    for i = 1:size(block, 1)
        % Extract the index for the current row
        index = indices(i);
        
        % Find the corresponding index in the pupil_combined data
        row_index = find(pupil_combined(:, 1) == index, 1);
        
        % Extract pupillometry data for the current row
        pupillometry_data{i} = pupil_combined(row_index, 2);
    end
end

% Average pupillometry data across all blocks
function avg_data = averageAcrossBlocks(data_blocks)
    % Find the maximum length of all blocks
    max_length = max(cellfun(@length, data_blocks));
    
    % Initialize a matrix to store interpolated data
    interpolated_data = nan(max_length, length(data_blocks));
    
    % Interpolate each block to the maximum length
    for i = 1:length(data_blocks)
        block_data = cell2mat(data_blocks{i});
        original_length = length(block_data);
        if original_length < max_length
            % Linear interpolation to the max length
            interpolated_data(:, i) = interp1(1:original_length, block_data, linspace(1, original_length, max_length));
        else
            % If block data is already at max length, just copy it
            interpolated_data(:, i) = block_data;
        end
    end
    
    % Calculate the mean across all blocks for each time point
    avg_data = mean(interpolated_data, 2, 'omitnan');
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
