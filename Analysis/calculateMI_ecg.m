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

    % Calculate MI for each time point within the window
    for t = 1:ntime_info
        if t <= window_size
            series = cdata(1:t);  % For initial points, use available data
        else
            series = cdata(t-window_size+1:t);  % Moving window
        end

        try
            % Calculate MI using mi_gg function
            MI(t) = mi_gg(series, series, false); % Placeholder for actual function
        catch ME
            warning('MI calculation failed at time point %d: %s', t, ME.message);
            MI(t) = NaN;
        end
    end
end
