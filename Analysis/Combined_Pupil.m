%% Load in eeglab
eeglab

%% Load in the data
% Load pupil data
pupil_data = load("//Client/D$/UGent_gerelateerd/Masterproef/Data/Pupilprepro/sub_032.mat");

% Load behavior file for pupil
beh_pupil = readtable("//Client/D$/UGent_gerelateerd/Masterproef/Data/RawData/sub-032/pupil/sub-032_task-memory_events.tsv", 'Delimiter', '\t', 'FileType', 'text');

%% Extract pupil diameter data
time = pupil_data.variables_to_export.world_index;
diameter = pupil_data.variables_to_export.diameter;

% Combine time and diameter into a single matrix
pupil_combined = [time, diameter];

%% Select and group trials (simplified without labels)
blocks_of_five = selectAndGroupTrials(beh_pupil, '5', 5);
blocks_of_nine = selectAndGroupTrials(beh_pupil, '9', 9);
blocks_of_thirteen = selectAndGroupTrials(beh_pupil, '3', 13);

%% Fill up missing indices in the blocks of trials
filled_blocks_of_five = cellfun(@fillMissingIndices, blocks_of_five, 'UniformOutput', false);
filled_blocks_of_nine = cellfun(@fillMissingIndices, blocks_of_nine, 'UniformOutput', false);
filled_blocks_of_thirteen = cellfun(@fillMissingIndices, blocks_of_thirteen, 'UniformOutput', false);


%% Extract pupillometry data for each condition and block
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
average_pupillometry_data = nanmean(cat(3, padded_data{:}), 3);

%% MI
