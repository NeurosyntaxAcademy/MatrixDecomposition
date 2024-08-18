% https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida

% CSV file containing the gene expression levels of 54676 genes (columns) from 151 samples (rows).
% There are 5 different types of breast cancer (plus healthy tissue) represented in this dataset 
% (column "type"). More information about this dataset, as well as other file formats such as 
% TAB and ARFF, data visualization, and classification and clustering benchmarks are 
% freely available at the official 
% CuMiDa website under the id GSE45827: http://sbcb.inf.ufrgs.br/cumida

% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy


%%
% change the directory
data = readtable(fullfile('C:\Users\HoffmanLab\Documents\Neurosyntax Academy\MatrixDecomposition\Principal Component Analysis', 'Breast_GSE45827.csv'));

Expression = table2array(data(:, 3:end));
Expression(isnan(Expression)) = 0;
CancerType = categorical(data.type);
CancerType = grp2idx(CancerType);

%%
% Standardize the data (zero mean, unit variance)
Expression = zscore(Expression, [], 2);

%%
% idx = randperm(size(Expression, 1));
% Expression = Expression(idx, :);

figure
imagesc(corr(Expression'))
colormap(perpl_RedBlueColormap), colorbar
caxis([-.5 .5])
pbaspect([1 1 1])
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)

% 41, 71, 72, 86, 93, 122

%%
% Perform PCA using the pca function
[coeff, score, latent, ~, explained] = pca(Expression);

% Plot the explained variance
figure;
plot(explained, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 30, 'LineWidth', 3);
title('Cumulative Explained Variance by Principal Components');
xlabel('Number of Principal Components');
ylabel('Explained Variance (%)');


%% Fit a robust regression line and determine the number of components
X = 1:length(explained);
b = robustfit(X,explained);
y_predicted = b(1) + b(2)*X;

figure;
plot(X, explained, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 30, 'LineWidth', 3); hold on
plot(X, y_predicted, 'Color', [1 0 0], 'Marker', '.', 'MarkerSize', 30, 'LineWidth', 3);
title('Cumulative Explained Variance by Principal Components');
xlabel('Number of Principal Components');
ylabel('Explained Variance (%)');


%%
% Select the number of principal components to retain (e.g., 2 for visualization)
num_pcs = 3;
reduced_data = score(:, 1:num_pcs);


% Plot the clustered brain areas in the reduced PCA space
figure;
gscatter(reduced_data(:, 1), reduced_data(:, 3), CancerType, [], [], 50);
title(['Explained Variance: ', num2str(sum(explained(1:2)))]);
xlabel('PC1');
ylabel('PC2');

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


gscatter3(reduced_data, CancerType)
title(['Explained Variance: ', num2str(sum(explained(1:3)))]);
xlabel('PC1');
ylabel('PC2');
zlabel('PC2');

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)



%% Interpreting loadings

imagesc(coeff(:, 1:3))
colormap(perpl_RedBlueColormap)

ylabel('Gene #'), xlabel('Components')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Sparse PCA
%[B SD L D paths] = spca(X, Gram, K, delta, stop, maxSteps, convergenceCriterion, verbose)

[B SD] = spca(Expression, [], 3, inf, -10000, 10000, 1e-9, true);
reduced_data = Expression*B;

gscatter3(reduced_data, CancerType)

title('Breast cancer gene expression - CuMiDa');
xlabel('sPC1');
ylabel('sPC2');
zlabel('sPC3');

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


imagesc(B)
colormap(perpl_RedBlueColormap)

ylabel('Gene #'), xlabel('Components')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)



%%

[reduced_data, umap, clusterIdentifiers]=run_umap(Expression,...
    'min_dist', 0.5, 'n_neighbors', 10, 'n_components', 3, ...
    'n_epochs', 5000, ...
    'metric', 'euclidean');


gscatter3(reduced_data, CancerType)

title('Breast cancer gene expression - CuMiDa');
xlabel('UMAP1');
ylabel('UMAP2');
zlabel('UMAP3');

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)



