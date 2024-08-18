% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy
% Toolbox: https://www2.imm.dtu.dk/projects/spasm/


%%

% Create synthetic data set
n = 1500; p = 500;
t = linspace(0, 1, p);
pc1 = max(0, (t - 0.5)> 0);
pc2 = 0.8*exp(-(t - 0.5).^2/5e-3);
pc3 = 0.4*exp(-(t - 0.15).^2/1e-3) + 0.4*exp(-(t - 0.85).^2/1e-3);
X = [ones(n/3,1)*pc1 + randn(n/3,p); ones(n/3,1)*pc2 + ...
    randn(n/3,p); ones(n/3,1)*pc3 + randn(n/3,p)];


% PCA and SPCA
[U S V] = svd(X, 'econ');
S = diag(S).^2/p; % PCA variances

% [B SD L D paths] = spca(X, Gram, K, delta, stop, maxSteps, convergenceCriterion, verbose)
[B SD L D] = spca(X, [], 3, inf, -[250 100 50], 3000, 1e-3, true);

figure(1)
plot(t, [pc1; pc2; pc3], 'LineWidth', 3); axis([0 1 -1.2 1.2]);
title('Noiseless data');
figure(2);
plot(t, X);  axis([0 1 -6 6]);
title('Data + noise');
figure(3);
plot(t, S(1:3).*(V(:,1:3)'), 'LineWidth', 3);  axis([0 1 -1.2 1.2]);
ylim([-10 20])
title('PCA');
figure(4)
plot(t, sqrt(SD)*ones(1,p).*(B'));  axis([0 1 -1.2 1.2]);
title('SPCA');


%%

figure

subplot(121)
imagesc(V(:,1:3))
cmap = perpl_RedBlueColormap();
colormap(cmap)
caxis([-0.1 0.1])

title('PCA [Explained Variance: 16.72%]')
xlabel('Component');
ylabel('Variables');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20, 'Xtick', 1:3)


subplot(122)
imagesc(B)
cmap = perpl_RedBlueColormap();
colormap(cmap)
caxis([-0.1 0.1])

title('Sparse PCA [Explained Variance: 16.41%]')
xlabel('Component');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20, 'Xtick', 1:3)


%% Parameter Tuning
X = zscore(X);

[coeff,score,latent] = pca(X, 'numcomponent', 3);

Beta = zeros(500, 3);
for i = 1:3
    %[b steps] = larsen(X, y, delta, stop, Gram, storepath, verbose)
    tmp = larsen(X, score(:, i), 1, 0, [], true, false);
    Beta(:, i) = tmp(:, end);
end


figure

subplot(121)
imagesc(Beta)
cmap = perpl_RedBlueColormap();
colormap(cmap)
caxis([-0.1 0.1])

title('Larse')
xlabel('Component');
ylabel('Variables');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20, 'Xtick', 1:3)


subplot(122)
imagesc(B)
cmap = perpl_RedBlueColormap();
colormap(cmap)
caxis([-0.1 0.1])

title('Sparse PCA [Explained Variance: 16.41%]')
xlabel('Component');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20, 'Xtick', 1:3)


%%


% Define the number of principal components to retain
K = 3;

% Define the lower and upper bounds for stop values
lb = [1 1 1]; % Lower bounds: 0
ub = [p p/2 p/4]; % Upper bounds: p (number of variables)

% Weights for combining variance difference and sparsity
weight_variance = 1; % Weight for variance difference
weight_sparsity = 0; % Weight for sparsity

% Objective function to minimize
objective = @(stop) SparsityVariance_objective(X, stop, K, weight_variance, weight_sparsity);


% Use the Genetic Algorithm to find the optimal stop values
options = optimoptions('ga', 'Display', 'iter', 'UseParallel', true);
[x,fval,exitflag,output,population,scores] = ga(objective, K, [], [], [], [], lb, ub, [], 1:K, options);

% Display the optimal stop values
disp('Optimal stop values:');
disp(optimal_stop);

scores_normalized = (scores - min(scores))/(max(scores) - min(scores));
figure,
imagesc(scores_normalized)
caxis([0 0.01])


