function MI = calculateMI_pupil(data)
    % Preprocess data to remove NaNs and Infs
    data = data(~isnan(data) & ~isinf(data));

    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = length(cdata);  % Assuming data is a single column

    % Initialize MI array
    MI = zeros(1, ntime_info);

    % Window size for calculating MI (you can adjust this size)
    window_size = 10;

    % Calculate MI for each time point within the window
    for t = 1:ntime_info
        if t <= window_size
            series = cdata(1:t);  % For initial points, use available data
        else
            series = cdata(t-window_size+1:t);  % Moving window
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

