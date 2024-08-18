classdef Preprocess
    properties
        random_seed
    end
    
    methods
        function obj = Preprocess(random_seed)
            if nargin < 1
                random_seed = [];
            end
            obj.random_seed = random_seed;
            if ~isempty(random_seed)
                rng(random_seed);
            end
        end
        
        function Z_corr = moe_correct_ridge(~, Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb)
            Z_corr = Z_orig;
            for i = 1:K
                Phi_Rk = Phi_moe .* R(i, :);
                x = Phi_Rk * Phi_moe' + lamb;
                W = (x \ Phi_Rk) * Z_orig';
                W(1, :) = 0; % do not remove the intercept
                Z_corr = Z_corr - (W' * Phi_Rk);
            end
            Z_cos = Z_corr / norm(Z_corr, 2);
        end
        
        function stdscale_quantile_celing(~, _adata, max_value, quantile_thresh)
            if nargin < 3
                max_value = [];
            end
            if nargin < 4
                quantile_thresh = [];
            end
            
            % Scaling
            _adata = scale_data(_adata, max_value);
            
            if ~isempty(quantile_thresh)
                if issparse(_adata.X)
                    threshval = quantile(full(_adata.X(:)), quantile_thresh);
                else
                    threshval = quantile(_adata.X(:), quantile_thresh);
                end
                _adata.X(_adata.X > threshval) = threshval;
            end
        end
        
        function make_count_hist(~, adata, num_cells)
            if nargin < 3
                num_cells = 1000;
            end
            
            z = full(adata.X(1:num_cells, :));
            y = z(:);
            histogram(y(y > 0), 100);
            title('Quantile thresholded normalized count distribution');
        end
        
        function _adata = filter_adata(~, _adata, filter_mito_thresh, min_cells_per_gene, min_counts_per_cell, filter_mito_genes, filter_dot_genes, makeplots)
            if nargin < 3
                filter_mito_thresh = [];
            end
            if nargin < 4
                min_cells_per_gene = 10;
            end
            if nargin < 5
                min_counts_per_cell = 500;
            end
            if nargin < 6
                filter_mito_genes = false;
            end
            if nargin < 7
                filter_dot_genes = true;
            end
            if nargin < 8
                makeplots = true;
            end
            
            if ~isempty(min_cells_per_gene)
                _adata = filter_genes(_adata, min_cells_per_gene);
            end
            
            _adata.obs.n_counts = sum(_adata.X, 2);
            
            if makeplots
                figure;
                histogram(log10(_adata.obs.n_counts), 100);
                title('log10 n_counts');
                ylim = get(gca, 'ylim');
                hold on;
                vline(log10(min_cells_per_gene), 'r');
                set(gca, 'ylim', ylim);
            end
            
            if ~isempty(min_counts_per_cell)
                _adata = filter_cells(_adata, min_counts_per_cell);
            end
            
            mt_genes = find(contains(_adata.var.index, 'MT-'));
            if ~isempty(filter_mito_thresh)
                num_mito = sum(_adata(:, mt_genes).X, 2);
                pct_mito = num_mito ./ _adata.obs.n_counts;
                _adata.obs.pct_mito = pct_mito;
                
                if makeplots
                    figure;
                    histogram(_adata.obs.pct_mito, 100);
                    title('pct_mito');
                end
                
                _adata = _adata(_adata.obs.pct_mito < filter_mito_thresh, :);
            end
            
            tofilter = [];
            if filter_dot_genes
                dot_genes = find(contains(_adata.var.index, '.'));
                tofilter = dot_genes;
            end
            
            if filter_mito_genes
                tofilter = [tofilter; mt_genes];
            end
            
            ind = ~ismember(1:length(_adata.var.index), tofilter);
            _adata = _adata(:, ind);
        end
        
        function [adata_RNA, tp10k, hvgs] = preprocess_for_cnmf(obj, _adata, feature_type_col, adt_feature_name, harmony_vars, n_top_rna_genes, librarysize_targetsum, max_scaled_thresh, quantile_thresh, makeplots, theta, save_output_base, max_iter_harmony)
            if nargin < 3
                feature_type_col = [];
            end
            if nargin < 4
                adt_feature_name = 'Antibody Capture';
            end
            if nargin < 5
                harmony_vars = [];
            end
            if nargin < 6
                n_top_rna_genes = 2000;
            end
            if nargin < 7
                librarysize_targetsum = 1e4;
            end
            if nargin < 8
                max_scaled_thresh = [];
            end
            if nargin < 9
                quantile_thresh = 0.9999;
            end
            if nargin < 10
                makeplots = true;
            end
            if nargin < 11
                theta = 1;
            end
            if nargin < 12
                save_output_base = [];
            end
            if nargin < 13
                max_iter_harmony = 20;
            end
            
            if ~iscell(_adata) && ~isempty(feature_type_col)
                adata_ADT = _adata(:, _adata.var.(feature_type_col) == adt_feature_name);
                adata_RNA = _adata(:, _adata.var.(feature_type_col) ~= adt_feature_name);
            elseif ~iscell(_adata)
                adata_RNA = _adata;
                adata_RNA.var.features_renamed = adata_RNA.var.index;
                adata_ADT = [];
            elseif length(_adata) == 2
                adata_RNA = _adata{1};
                adata_ADT = _adata{2};
                if size(adata_ADT, 1) ~= size(adata_RNA, 1)
                    error("ADT and RNA AnnDatas don't have the same number of cells");
                elseif any(adata_ADT.obs.index ~= adata_RNA.obs.index)
                    error("Inconsistency of the index for the ADT and RNA AnnDatas");
                end
            else
                error('data should either be an AnnData object or a list of 2 AnnData objects');
            end
            
            tp10k = normalize_per_cell(adata_RNA, librarysize_targetsum);
            [adata_RNA, hvgs] = obj.normalize_batchcorrect(adata_RNA, harmony_vars, n_top_rna_genes, librarysize_targetsum, max_scaled_thresh, quantile_thresh, theta, makeplots, max_iter_harmony);
            
            if ~isempty(adata_ADT)
                adata_ADT = adata_ADT(adata_RNA.obs.index, :);
                adata_ADT = normalize_per_cell(adata_ADT, librarysize_targetsum);
                
                merge_var = [tp10k.var; adata_ADT.var];
                tp10k = struct('X', [tp10k.X, adata_ADT.X], 'obs', tp10k.obs, 'var', merge_var);
            end
            
            if ~isempty(save_output_base)
                save([save_output_base '.Corrected.HVG.Varnorm.mat'], 'adata_RNA');
                save([save_output_base '.TP10K.mat'], 'tp10k');
                writematrix(hvgs, [save_output_base '.Corrected.HVGs.txt']);
            end
        end
        
        function [adata, hvgs] = normalize_batchcorrect(obj, _adata, harmony_vars, n_top_genes, librarysize_targetsum, max_scaled_thresh, quantile_thresh, theta, makeplots, max_iter_harmony)
            if nargin < 3
                harmony_vars = [];
            end
            if nargin < 4
                n_top_genes = 2000;
            end
            if nargin < 5
                librarysize_targetsum = 1e4;
            end
            if nargin < 6
                max_scaled_thresh = [];
            end
            if nargin < 7
                quantile_thresh = 0.9999;
            end
            if nargin < 8
                theta = 1;
            end
            if nargin < 9
                makeplots = true;
            end
            if nargin < 10
                max_iter_harmony = 20;
            end
            
            if ~isempty(n_top_genes)
                _adata = highly_variable_genes(_adata, 'seurat_v3', n_top_genes);
            elseif ~ismember('highly_variable', _adata.var.Properties.VariableNames)
                error("If a numeric value for n_top_genes is not provided, you must include a highly_variable column in _adata");
            end
            
            if ~isempty(harmony_vars)
                anorm = normalize_per_cell(_adata, librarysize_targetsum);
                anorm = anorm(:, _adata.var.highly_variable);
                obj.stdscale_quantile_celing(anorm, max_scaled_thresh, quantile_thresh);
                
                _adata = _adata(:, _adata.var.highly_variable);
                obj.stdscale_quantile_celing(_adata, max_scaled_thresh, quantile_thresh);
                
                if makeplots
                    obj.make_count_hist(anorm, 1000);
                end
                
                _adata = pca(_adata, true);
                if makeplots
                    pca_variance_ratio(_adata, true, 50);
                end
                
                _adata.obsm.X_pca = anorm.obsm.X_pca;
                
                if ~normalize_librarysize
                    [~, _adata.obsm.X_pca_harmony] = obj.harmony_correct_X(full(_adata.X), _adata.obs, _adata.obsm.X_pca, harmony_vars, theta, max_iter_harmony);
                else
                    [~, _adata.obsm.X_pca_harmony] = obj.harmony_correct_X(full(anorm.X), anorm.obs, anorm.obsm.X_pca, harmony_vars, theta, max_iter_harmony);
                end
            else
                if normalize_librarysize
                    _adata = normalize_per_cell(_adata, librarysize_targetsum);
                end
                
                _adata = _adata(:, _adata.var.highly_variable);
                obj.stdscale_quantile_celing(_adata, max_scaled_thresh, quantile_thresh);
                if makeplots
                    obj.make_count_hist(_adata, 1000);
                end
            end
            
            hvgs = _adata.var.index;
        end
        
        function [X_corr, X_pca_harmony] = harmony_correct_X(~, X, obs, pca, harmony_vars, theta, max_iter_harmony)
            try
                harmonypy;
            catch
                error("harmonypy is not installed. Please install it using 'pip install harmonypy' before proceeding.");
            end
            
            harmony_res = harmonypy.run_harmony(pca, obs, harmony_vars, 'max.iter.harmony', max_iter_harmony, 'theta', theta);
            
            X_pca_harmony = harmony_res.Z_corr';
            [~, X_corr, ~, ~] = moe_correct_ridge(X', [], [], harmony_res.R, [], harmony_res.K, [], harmony_res.Phi_moe, harmony_res.lamb);
            X_corr = X_corr';
            
            X_corr(X_corr < 0) = 0;
        end
        
        function _adata = select_features_MI(obj, _adata, cluster, max_scaled_thresh, quantile_thresh, n_top_features, makeplots)
            if nargin < 4
                max_scaled_thresh = [];
            end
            if nargin < 5
                quantile_thresh = 0.9999;
            end
            if nargin < 6
                n_top_features = 70;
            end
            if nargin < 7
                makeplots = true;
            end
            
            _adata = normalize_per_cell(_adata);
            obj.stdscale_quantile_celing(_adata, max_scaled_thresh, quantile_thresh);
            
            if issparse(_adata.X)
                res = mutual_info_classif(full(_adata.X), cluster);
            else
                res = mutual_info_classif(_adata.X, cluster);
            end
            
            res = sortrows(table(res, _adata.var.index), 'res', 'descend');
            resdf = [res(:, 2), res(:, 1), diff(res(:, 1))];
            
            if makeplots
                figure;
                scatter(1:length(resdf), resdf(:, 2));
                ylabel('MI');
                xlabel('MI Rank');
                ylim = get(gca, 'ylim');
                hold on;
                vline(n_top_features, 'k--');
                set(gca, 'ylim', ylim);
            end
            
            _adata.var.MI = resdf(:, 2);
            _adata.var.MI_Rank = (1:length(resdf))';
            _adata.var.MI_diff = [NaN; diff(_adata.var.MI)];
            _adata.var.highly_variable = _adata.var.MI_Rank < n_top_features;
        end
    end
end

function data = scale_data(data, max_value)
    % Add your custom scaling function here
    if ~isempty(max_value)
        data.X = data.X / max(max(data.X)) * max_value;
    end
end

function data = normalize_per_cell(data, librarysize_targetsum)
    % Add your custom normalization function here
    data.X = data.X / sum(data.X, 2) * librarysize_targetsum;
end

function data = filter_genes(data, min_cells_per_gene)
    % Add your custom gene filtering function here
    gene_counts = sum(data.X > 0, 1);
    data = data(:, gene_counts >= min_cells_per_gene);
end

function data = filter_cells(data, min_counts_per_cell)
    % Add your custom cell filtering function here
    cell_counts = sum(data.X, 2);
    data = data(cell_counts >= min_counts_per_cell, :);
end

function data = highly_variable_genes(data, flavor, n_top_genes)
    % Add your custom highly variable gene selection function here
    % This is a placeholder implementation
    gene_variance = var(data.X, 0, 1);
    [~, sorted_idx] = sort(gene_variance, 'descend');
    data.var.highly_variable = false(height(data.var), 1);
    data.var.highly_variable(sorted_idx(1:n_top_genes)) = true;
end

function pca(data, use_highly_variable)
    % Add your custom PCA function here
    % This is a placeholder implementation
    if use_highly_variable
        data.X = data.X(:, data.var.highly_variable);
    end
    [coeff, score] = pca(data.X);
    data.obsm.X_pca = score;
end

function pca_variance_ratio(data, log, n_pcs)
    % Add your custom PCA variance ratio plot function here
    % This is a placeholder implementation
    pca_var = var(data.obsm.X_pca, 0, 1);
    pca_var = pca_var / sum(pca_var);
    if log
        semilogy(pca_var(1:n_pcs));
    else
        plot(pca_var(1:n_pcs));
    end
    title('PCA variance ratio');
end

function mutual_info_classif(X, y)
    % Add your custom mutual information classification function here
    % This is a placeholder implementation
    res = mutualinfo(X, y);
end

function vline(x, line_style)
    % Add your custom vertical line plotting function here
    % This is a placeholder implementation
    line([x x], ylim, 'LineStyle', line_style, 'Color', 'k');
end
