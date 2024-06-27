function II = calculateII_different(MI, JMI)
    % Calculate Interaction Information (II) based on MI and JMI
    ntime_info = length(MI); % number of time points
    II = NaN(ntime_info, ntime_info); % Initialize II matrix with NaNs
    
    % Compute II directly using JMI vector
    idx = 1; % Start index for JMI vector
    for t1 = 1:ntime_info
        for t2 = 1:ntime_info
            % Check for NaNs in MI or JMI values
            if isnan(MI(t1)) || isnan(MI(t2)) || idx > length(JMI)
                II(t1, t2) = NaN;
            else
                % Calculate II using the formula: II(t1, t2) = JMI(t1, t2) - MI(t1) - MI(t2)
                II(t1, t2) = JMI(idx) - MI(t1) - MI(t2);
                idx = idx + 1; % Move to the next index in JMI vector
            end
        end
    end
end

