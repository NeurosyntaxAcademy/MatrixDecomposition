%%% Example: Happiness, trust, and deaths under COVID-19
% https://worldhappiness.report/ed/2021/happiness-trust-and-deaths-under-covid-19/
% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy


%%
originalTable = readtable('world-happiness-report-2021.csv');
Data = originalTable(:, [8 9 10 11 12]);
Data = table2array(Data);

%%

% Assume Data is the original data matrix where rows are observations and columns are features

% Standardize the data (each column has mean 0 and variance 1)
standardizedData = zscore(Data);

% Compute the covariance matrix
covMatrix = corr(standardizedData);

% Perform eigen decomposition of the covariance matrix
[eigenvectors, eigenvalues] = eig(covMatrix);

% Extract the eigenvalues from the diagonal matrix
eigenvalues = diag(eigenvalues);

% Sort the eigenvalues and corresponding eigenvectors in descending order
[eigenvalues, order] = sort(eigenvalues, 'descend');
eigenvectors = eigenvectors(:, order);

%% Explained Variance
figure('units','normalized','outerposition',[0 0 1 1])

plot(eigenvalues/sum(eigenvalues), 'LineWidth', 3, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 100)

title('Scree Plot - Explained Variance')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%%
figure('units','normalized','outerposition',[0 0 1 1])

stem(eigenvectors(:, 2), 'Color', [0 0 0], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 50)
title('World Happiness Scores')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
set(gca, 'XTick', [1:5], 'XTicklabel', {'Social', 'Life', 'Choice', 'Generosity', 'Corruption'})

xlim([0 6])




