function MI = calculateMI_pupil(data)
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

        try
            % Center the data
            centered_series = series - mean(series); 
            
            % Calculate covariance matrix
            cov_matrix = cov(centered_series);

            if isscalar(cov_matrix)
                eigenvalues = cov_matrix; % For scalar, eigenvalue is the value itself
            else
                eigenvalues = eig(cov_matrix);
            end

            % Check if the covariance matrix is positive definite
            if all(eigenvalues > 0)
                % Handle scalar case separately
                if isscalar(cov_matrix)
                    mi_value = -0.5 * log(2 * pi * exp(1) * cov_matrix); % Simplified MI for Gaussian
                else
                    % Calculate MI using mi_gg function for non-scalar cases
                    mi_value = mi_gg(centered_series, centered_series, false); % Placeholder for actual function
                end
                MI(t) = mi_value;
            else
                warning('Covariance matrix is not positive definite at time point %d.', t);
                MI(t) = NaN;
            end
        catch ME
            warning('MI calculation failed at time point %d: %s', t, ME.message);
            MI(t) = NaN;
        end
    end
end