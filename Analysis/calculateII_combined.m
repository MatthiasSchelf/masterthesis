function II = calculateII_combined(MI, JMI)
    % Calculate Interaction Information (II) based on MI and JMI
    ntime_info = length(MI); % number of time points
    II = NaN(ntime_info, ntime_info); % Initialize II matrix with NaNs
    
    % Compute II directly using MI and JMI vectors
    for t1 = 1:ntime_info
        for t2 = 1:ntime_info
            % Check for NaNs in MI or JMI values
            if isnan(MI(t1)) || isnan(MI(t2)) || isnan(JMI(t1)) || isnan(JMI(t2))
                II(t1, t2) = NaN;
            else
                % Calculate II using the formula: II(t1, t2) = JMI(t1) - MI(t1) - MI(t2)
                % Note: This assumes JMI is provided for each time point pair (t1, t2)
                II(t1, t2) = JMI(t1) - MI(t1) - MI(t2);
            end
        end
    end
end


