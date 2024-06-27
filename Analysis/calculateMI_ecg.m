function MI = calculateMI_ecg(data)
    % Preprocess data to remove NaNs and Infs
    data = data(~isnan(data) & ~isinf(data));

    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = length(cdata);  % Assuming data is a single column

    % Initialize MI array
    MI = zeros(1, ntime_info);

    % Window size for calculating MI (you can adjust this size)
    window_size = 10;

    % Regularization term
    regularization_term = 1e-6;

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

            % Check if the covariance matrix is positive definite
            [~, p] = chol(cov_matrix + regularization_term * eye(size(cov_matrix)));

            if p ~= 0
                % Covariance matrix is not positive definite
                % Apply regularization until it becomes positive semi-definite
                epsilon = 1e-10;
                while p ~= 0
                    cov_matrix = cov_matrix + epsilon * eye(size(cov_matrix));
                    [~, p] = chol(cov_matrix);
                    epsilon = epsilon * 10; % Increase epsilon if necessary
                end
            end

            % Calculate MI using mi_gg function
            if isscalar(cov_matrix)
                mi_value = -0.5 * log(2 * pi * exp(1) * cov_matrix); % Simplified MI for Gaussian
            else
                mi_value = mi_gg(centered_series, centered_series, false); % Placeholder for actual function
            end

            MI(t) = mi_value;
        catch ME
            warning('MI calculation failed at time point %d: %s', t, ME.message);
            MI(t) = NaN;
        end
    end
end
