function MI = calculateMI_ecg(data)
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

        % Calculate MI using mi_gg function
        try
            MI(t) = mi_gg(series, series, false); % Placeholder
        catch ME
            warning('MI calculation failed at time point %d: %s', t, ME.message);
            MI(t) = NaN;
        end
    end
end