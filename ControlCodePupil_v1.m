%%  Step 1: Load in data

data = readtable('sub-032_task-memory_pupil.tsv', 'Delimiter', '\t', 'FileType', 'text');

%% Step 2: Visualize raw data

% Example data (replace this with your actual variable names)
time = data.world_index;
pupil_size = data.diameter;

subset_time = time(1:10000);
subset_pupil_size = pupil_size(1:10000);

% Plotting the raw pupillometry data for the first 2000 time points
figure;
plot(subset_time, subset_pupil_size);
title('Pupillometry Data');
xlabel('Time (seconds)');
ylabel('Pupil Size');

%% Step 3: Identify and Remove Blinks using Z-Score

% Define z-score threshold for blink detection (adjust as needed)
z_threshold = 1.7; % You can experiment with different threshold values

% Calculate z-scores for the pupil size
z_scores = zscore(subset_pupil_size);

% Find indices of blinks based on z-score threshold
blink_indices = find(abs(z_scores) > z_threshold);

% Remove blinks from the data
cleaned_time = subset_time;
cleaned_pupil_size = subset_pupil_size;
cleaned_pupil_size(blink_indices) = NaN; % Set blink values to NaN or remove them as needed

% Plot the cleaned data with y-axis limits set to 0-60
figure;
plot(cleaned_time, cleaned_pupil_size);
title('Cleaned Pupillometry Data (Blinks Removed)');
xlabel('Time (seconds)');
ylabel('Pupil Size');

% Set y-axis limits to 0-60
ylim([0 60]);

%% Step 4: Interpolate the Removed Blink Data

% Find NaN values (indicating removed blinks)
nan_indices = isnan(cleaned_pupil_size);

% Interpolate the removed blink data using fillmissing
interpolated_pupil_size = fillmissing(cleaned_pupil_size, 'linear');

% Plot the data with interpolated blinks
figure;
plot(cleaned_time, interpolated_pupil_size);
title('Interpolated Pupillometry Data');
xlabel('Time (seconds)');
ylabel('Pupil Size');

% Optionally, you can replace the NaN values in the original cleaned data
% with the interpolated values
cleaned_pupil_size(nan_indices) = interpolated_pupil_size(nan_indices);

% Plot the final cleaned and interpolated data
figure;
plot(cleaned_time, cleaned_pupil_size);
title('Final Cleaned and Interpolated Pupillometry Data');
xlabel('Time (seconds)');
ylabel('Pupil Size');
% Set y-axis limits to 0-60
ylim([0 60]);

%% Step 5: Baseline Correction

% Calculate baseline using mean or median of the cleaned and interpolated data
baseline = mean(interpolated_pupil_size, 'omitnan'); % You can also use median()

% Subtract the baseline from the interpolated pupil size
baseline_corrected_pupil_size = interpolated_pupil_size - baseline;

% Plot the baseline-corrected data
figure;
plot(cleaned_time, baseline_corrected_pupil_size);
title('Baseline-Corrected Pupillometry Data');
xlabel('Time (seconds)');
ylabel('Baseline-Corrected Pupil Size');









