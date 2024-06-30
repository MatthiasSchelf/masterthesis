% Function to get the difference of pupil data between the multivariate
% baseline and periods of interest
function multivariate_difference_pupil = multivariate_pupil_difference(multivariate_time_series_pupil, single_baseline_value_pupil)
    % Extract mu values
    multivariate_mu_pupil = cellfun(@(x) x.mu, multivariate_time_series_pupil, 'UniformOutput', false);
    
    % Subtract the one baseline value from each separate time point in the period of interest 
    multivariate_difference_pupil_temp = cellfun(@(mu) mu - single_baseline_value_pupil, multivariate_mu_pupil, 'UniformOutput', false);

    % Convert the cell array to a double array and concatenate all
    % differences into a single column (just for further calculation
    % purposes)
    multivariate_difference_pupil = [];
    for i = 1:length(multivariate_difference_pupil_temp)
        mu_diff = multivariate_difference_pupil_temp{i};
        % Reshape mu_diff to a column vector and concatenate
        multivariate_difference_pupil = [multivariate_difference_pupil; mu_diff(:)];
    end
end


