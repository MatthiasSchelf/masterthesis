% Function that calculates the difference between the baseline and perios
% of interest for EEG
function multivariate_difference_eeg = multivariate_eeg_difference(multivariate_time_series_eeg, single_baseline_value_eeg)
    % Extract mu values
    multivariate_mu_pupil = cellfun(@(x) x.mu, multivariate_time_series_eeg, 'UniformOutput', false);
    
    % Subtract the single baseline value from each timepoint of the period
    % of interest
    multivariate_difference_pupil = cellfun(@(mu) mu - single_baseline_value_eeg, multivariate_mu_pupil, 'UniformOutput', false);

    % Ensure all mu values are vectors and concatenate them into a single column
    multivariate_difference_eeg = [];
    for i = 1:length(multivariate_difference_pupil)
        mu_diff = multivariate_difference_pupil{i};
        % Reshape mu_diff to a column vector and concatenate (easier for further
        % calculations).
        multivariate_difference_eeg = [multivariate_difference_eeg; mu_diff(:)];
    end
end
