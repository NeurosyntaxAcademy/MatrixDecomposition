%% https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_pbmc_example_data.ipynb
% Data: http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

%% file
clear all
cd('C:\Users\abbasps\Documents\Matrix Decomposition\pbmc3k_filtered_gene_bc_matrices\filtered_gene_bc_matrices\hg19')
filename = 'matrix.mtx';

%% read Matrix
[A,rows,cols,entries,rep,field,symm] = mmread('matrix.mtx');
% genes = tdfread('genes.tsv');
% barcodes = tdfread('barcodes.tsv');
A = full(A);

% [A, rowNames, colNames] = readMtxFile('matrix.mtx');

%% Preprocessing
% For each dataset, we removed cells with fewer than 1000 unique molecular identifiers (UMIs)
% detected. We also filtered out genes that were not detected in at least 1 out of 500 cells.


% filter cells with fewer than 200 genes
Unique_UMI_per_cell = sum(logical(A), 1);  % Sum of counts per cell (column-wise sum)
cells_to_keep = Unique_UMI_per_cell >= 200;
A_filtered_cells = A(:, cells_to_keep);

% This is a weaker threshold than above. It is just to population the n_counts column in adata
UMI_counts_per_cell = sum(A_filtered_cells, 1);  % Sum of counts per cell (column-wise sum)
cells_to_keep = UMI_counts_per_cell >= 200;
A_filtered_cells = A_filtered_cells(:, cells_to_keep);

% filter genes detected in fewer than 3 cells
gene_detection_count = sum(logical(A_filtered_cells) > 0, 2);  % Number of cells in which each gene is detected
genes_to_keep = gene_detection_count >= 3;
A_filtered = A_filtered_cells(genes_to_keep, :);

% Display the size of the filtered matrix
disp(['Original matrix size: ', num2str(size(A))]);
disp(['Filtered matrix size: ', num2str(size(A_filtered))]);

%%

histogram(log10(sum(A_filtered)), 100, 'FaceColor', [0 0 0], 'FaceAlpha', 1)
xlabel('log10 Counts Per Cell');
ylabel('# Cells');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

%% V-Score
% Klein AMMazutis LAkartuna ITallapragada NVeres ALi VPeshkin LWeitz DAKirschner MW (2015)
% Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells Cell
% https://doi.org/10.1016/j.cell.2015.04.044

nu = vscore(A_filtered);
H = 2000;  % Number of top genes to select
[~, top_gene_indices] = maxk(nu, H);


% Filter the matrix to keep only the top H genes
A_top_genes = A_filtered(top_gene_indices, :);

% Display the size of the matrix with top genes
disp(['Matrix size with top ', num2str(H), ' genes: ', num2str(size(A_top_genes))]);

% The matrix A_top_genes now contains only the top H most over-dispersed genes

A_top_genes = zscore(A_top_genes, [], 2);
A_top_genes = A_top_genes';

size(A_top_genes)

% [From Paper] Note that we do not perform any cell count normalization
% This is because cells with more counts can contribute more information to the model.

%%
A_top_genes = full(A_top_genes);

[U,S,V] = svd(A_top_genes, 'econ');
S = diag(S);
score = A_top_genes * V(:, 1:15);

figure;
plot(S, 'ko-', 'LineWidth', 2, 'MarkerSize', 10);
title('Singular Values');
xlabel('Component');
ylabel('Value');
grid on;
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)


%% UMAP for visualization

% [reduction, umap, clusterIdentifiers, extras]=run_umap()
[Ared, umap, clusterIdentifiers]=run_umap(score,...
    'min_dist', 0.5, 'n_neighbors', 50, 'n_components', 2, ...
    'n_epochs', 1000, ...
    'metric', 'euclidean');

scatter(Ared(:, 1), Ared(:, 2), 50, [0 0 0], 'MarkerFaceAlpha', .5)
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

%%

Options = [];
Options = statset('MaxIter', 10000, 'TolFun', 1e-8, 'UseParallel', false, 'TolX', 1e-8);

numComponents = 15; % Set the number of components for the factorization
[W, H] = nnmf(A_top_genes, numComponents, 'options', Options);


%%

% [reduction, umap, clusterIdentifiers, extras]=run_umap()
[Wred, umap, clusterIdentifiers]=run_umap(W,...
    'min_dist', 0.5, 'n_neighbors', 50, 'n_components', 2, ...
    'n_epochs', 1000, ...
    'metric', 'euclidean');

scatter(Wred(:, 1), Wred(:, 2), 50, [0 0 0], 'filled', 'MarkerFaceAlpha', 0.5)
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

%%
% Colormaps from: https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap

imagesc(H(1:5, :)')
cmap = slanCM('amethyst');
colormap(cmap), colormap
caxis([0 .1])

xlabel('Components');
ylabel('Gene');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)






