function single_baseline_value = multivariate_reduce_one(time_resolved_distributions)
    % Initialize the sum of mean values and a counter for the number of mean values
    sum_of_mu_values = 0;
    num_mu_values = 0;

    % Iterate through each time point's distribution
    for t = 1:length(time_resolved_distributions)
        % Extract the mean value (mu)
        mu = time_resolved_distributions{t}.mu;

        % Sum the mean values (assuming mu is a row vector)
        sum_of_mu_values = sum_of_mu_values + sum(mu);
        
        % Count the number of mean values
        num_mu_values = num_mu_values + length(mu);
    end

    % Compute the average of all mean values
    single_baseline_value = sum_of_mu_values / num_mu_values;
end

