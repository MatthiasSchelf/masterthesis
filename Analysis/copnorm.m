% Copula normalization function
function cn_data = copnorm(data)
    % Rank transform the data
    [~, ranks] = sort(data);
    [~, ranks] = sort(ranks);
    ranks = ranks / (length(data) + 1);
    
    % Apply the inverse normal (Gaussian) CDF
    cn_data = norminv(ranks);
end