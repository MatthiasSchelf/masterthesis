% Function to make the multivariate baseline one single number, so that it
% can be used to substract from the period of interest to create the change
% from baseline.

function single_baseline_value = multivariate_reduce_one(time_resolved_distributions)
    
    % Initialize the sum of multivariate values and a counter for the number of multivariate values
    sum_of_mu_values = 0;
    num_mu_values = 0;

    % Iterate through each time point's distribution
    for t = 1:length(time_resolved_distributions)
        % Extract the mean value (mu)
        mu = time_resolved_distributions{t}.mu;

        % Sum the mean values
        sum_of_mu_values = sum_of_mu_values + sum(mu);
        
        % Count the number of mean values
        num_mu_values = num_mu_values + length(mu);
    end

    % Average of multivatiate baseline in one number.
    single_baseline_value = sum_of_mu_values / num_mu_values;
end

