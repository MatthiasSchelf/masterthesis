function multivariate_difference_ecg = multivariate_ecg_difference(multivariate_time_series_ecg, single_baseline_value_ecg)
    % Extract mu values
    multivariate_mu_ecg = cellfun(@(x) x.mu, multivariate_time_series_ecg, 'UniformOutput', false);
    
    % Subtract single baseline value from each mu
    multivariate_difference_ecg_temp = cellfun(@(mu) mu - single_baseline_value_ecg, multivariate_mu_ecg, 'UniformOutput', false);

    % Convert the cell array to a double array and concatenate all differences into a single column
    multivariate_difference_ecg = [];
    for i = 1:length(multivariate_difference_ecg_temp)
        mu_diff = multivariate_difference_ecg_temp{i};
        % Reshape mu_diff to a column vector and concatenate
        multivariate_difference_ecg = [multivariate_difference_ecg; mu_diff(:)];
    end
end
