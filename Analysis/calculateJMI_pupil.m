function JMI = calculateJMI_pupil(data)
    % Preprocess data to remove NaNs and Infs
    data = data(~isnan(data) & ~isinf(data));

    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = length(cdata);  % Assuming data is a single column

    % Initialize JMI array
    JMI = zeros(1, ntime_info);

    % Window size for calculating JMI (you can adjust this size)
    window_size = 10;

    % Regularization term
    regularization_term = 1e-6;

    % Calculate JMI for each time point within the window
    for t = 1:ntime_info
        if t <= window_size
            series = cdata(1:t);  % For initial points, use available data
        else
            series = cdata(t-window_size+1:t);  % Moving window
        end

        try
            % Calculate variance with regularization
            var_series = var(series) + regularization_term;

            % Check if the variance is positive
            if var_series > 0
                % Calculate JMI using mi_gg function for Gaussian cases
                % Note: mi_gg is a placeholder for actual function
                mi_value = mi_gg(series, series, false);
                JMI(t) = mi_value;
            else
                warning('Variance is not positive at time point %d.', t);
                JMI(t) = NaN;
            end
        catch ME
            warning('JMI calculation failed at time point %d: %s', t, ME.message);
            JMI(t) = NaN;
        end
    end
end

