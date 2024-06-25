function copnormed_data = copnorm(data)
    % Copula normalization function
    n = length(data);
    ranks = tiedrank(data);
    copnormed_data = (ranks - 0.5) / n;
end