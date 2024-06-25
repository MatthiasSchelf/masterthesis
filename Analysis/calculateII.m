function II = calculateII(MI, JMI)
    % Calculate Interaction Information (II) based on MI and JMI
    ntime_info = size(MI, 2);
    II = zeros(ntime_info, ntime_info);
    
    % Loop through pairs of time points
    for t1 = 1:ntime_info
        for t2 = (t1 + 1):ntime_info
            % Check for NaNs or Infs in the MI values
            if isnan(MI(t1)) || isnan(MI(t2)) || isnan(JMI(t1, t2))
                II(t1, t2) = NaN;
                continue;
            end
            
            % Calculate II using the formula: II(t1, t2) = JMI(t1, t2) - MI(t1) - MI(t2)
            II(t1, t2) = JMI(t1, t2) - MI(t1) - MI(t2);
        end
    end
    
    % Symmetrize the II matrix
    II = II + II';
end