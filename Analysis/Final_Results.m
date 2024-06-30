% Load the saved participant results and run this code to combine all
% results of II together for every participant (this code has noe been
% tested).

% List of participant IDs
participants = {'sub-032'};  % Add more participant IDs as needed
participant_field_names = cellfun(@(x) strrep(x, '-', '_'), participants, 'UniformOutput', false);

% Initialize variables to store the sum of II matrices
sum_II_average = [];
sum_II_multivariate = [];

% Loop over each participant to accumulate the II matrices
for i = 1:length(participant_field_names)
    participant_field_name = participant_field_names{i};
    
    % Extract II matrices
    II_average = participant_results.(participant_field_name).II_average;
    II_multivariate = participant_results.(participant_field_name).II_multivariate;
    
    % Initialize sum matrices if empty
    if isempty(sum_II_average)
        sum_II_average = zeros(size(II_average));
    end
    if isempty(sum_II_multivariate)
        sum_II_multivariate = zeros(size(II_multivariate));
    end
    
    % Accumulate the II matrices
    sum_II_average = sum_II_average + II_average;
    sum_II_multivariate = sum_II_multivariate + II_multivariate;
end

% Calculate the average II matrices
num_participants = length(participant_field_names);
average_II_average = sum_II_average / num_participants;
average_II_multivariate = sum_II_multivariate / num_participants;

% Display the averaged II matrices
disp('Averaged II (Average):');
disp(average_II_average);

disp('Averaged II (Multivariate):');
disp(average_II_multivariate);

