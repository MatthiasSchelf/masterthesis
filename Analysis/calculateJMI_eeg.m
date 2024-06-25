function JMI = calculateJMI_eeg(data)
    % Transform data using copula normalization
    cdata = copnorm(data);
    ntime_info = size(cdata, 2);

    % Initialize JMI array
    JMI = zeros(ntime_info, ntime_info);

    % Calculate JMI for each pair of time points
    for t1 = 1:ntime_info
        for t2 = (t1 + 1):ntime_info
            % Extract the time series at time t1 and t2
            series1 = cdata(:, t1);
            series2 = cdata(:, t2);

            % Check for NaNs or Infs
            if any(isnan(series1)) || any(isinf(series1)) || any(isnan(series2)) || any(isinf(series2))
                warning('Data contains NaNs or Infs at time points %d and %d. Skipping.', t1, t2);
                JMI(t1, t2) = NaN;
                continue;
            end

            % Calculate JMI using mi_gg function
            try
                JMI(t1, t2) = mi_gg([series1 series2], [series1 series2], false); % Placeholder
            catch ME
                warning('JMI calculation failed at time points %d and %d: %s', t1, t2, ME.message);
                JMI(t1, t2) = NaN;
            end
        end
    end
end
