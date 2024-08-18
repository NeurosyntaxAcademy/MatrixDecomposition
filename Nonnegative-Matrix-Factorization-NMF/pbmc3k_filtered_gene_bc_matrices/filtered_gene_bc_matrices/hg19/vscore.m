function nu = vscore(data) %alpha
    % vscore Identifies highly variable genes from expression data
    %
    % Inputs:
    %   data - A matrix where rows represent genes and columns represent samples
    %   alpha - The significance level for FDR control (typically 0.05)
    %
    % Output:
    %   significant_genes - Indices of genes that are significantly variable

    % Calculate mean and CV for each gene
    m = mean(data, 2); % Row-wise mean
    sd = std(data, 0, 2); % Standard deviation for each gene, 0 for default normalization
    cv = sd ./ m; % Coefficient of variation

    % Calculate test statistic ν for each gene
    cvsq = cv.^2;
    nu = (cvsq ./ (1 + cvsq)) ./ (m + cvsq);

    % Define null model parameters for Poisson distribution, if known
    % lambda = m; % Assuming the mean of the data approximates the lambda of the Poisson distribution
    % cv_null = sqrt(1 ./ lambda); % CV under Poisson assumption

    % Calculate z-scores
    % z = (cv - cv_null) ./ std(cv_null);

    % % Calculate p-values
    % p = normcdf(-abs(z), 0, 1) * 2; % Two-tailed test
    % 
    % % Apply Benjamini and Hochberg FDR control
    % [sorted_p, idx] = sort(p);
    % adjusted_p = cumsum(sorted_p) ./ (1:numel(sorted_p))';
    % significant = sorted_p < adjusted_p * alpha;
    % 
    % % Find significant genes
    % significant_genes = find(significant(idx));

    % Return the indices of significant genes
end
