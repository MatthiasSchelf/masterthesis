% Load ECG data
load("sub_032.mat")
% Load pupil diameter data
load("sub-032_task-memory_ecg.mat")

% Extract pupil diameter data
pupil_data = variables_to_export;
time = pupil_data.world_index;
diameter = pupil_data.diameter;

% Combine time and diameter into a single matrix
pupil_combined = [time, diameter];

% Extract ECG data from the HEP structure
ecg_data = HEP.ecg;

% Check the number of trials (samples) in both ECG data and pupil data
num_trials_x = size(ecg_data, 1); % Assuming each row is a trial/sample
num_trials_y = size(pupil_combined, 1);

% If the number of trials doesn't match, adjust the data accordingly
if num_trials_x ~= num_trials_y
    % Adjust the data to match the number of trials
    min_trials = min(num_trials_x, num_trials_y);
    ecg_data = ecg_data(1:min_trials, :);
    pupil_combined = pupil_combined(1:min_trials, :);
end

% Now you can use the gcmi_cc function
gcmi_score = gcmi_cc(ecg_data, pupil_combined);
disp(['GCMi Score: ', num2str(gcmi_score)]);

