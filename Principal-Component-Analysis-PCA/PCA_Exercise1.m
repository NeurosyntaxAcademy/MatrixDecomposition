% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy

%%

clear all, clc

%% Simulate a dataset
Data = zeros(1000, 6);
Data(randi(1000, [1 100]), 1:3) = 1;

Data(randi(1000, [1 100]), 4:6) = 1;
Data(Data(:, 1) == 1, 4:6) = 0;


% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

imagesc(Data')
colormap(flip(gray(256)))
xlabel('Sample (e.g. Time)')
ylabel('Variable (e.g. Neuron/Gene)')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Compute correlation matrix

% Method 1
standardizedData = (Data - mean(Data))./std(Data);
corrMatrix = (standardizedData' * standardizedData) / (size(Data, 1) - 1);


% Method 2
% standardizedData = zscore(standardizedData);
% corrMatrix = cov(standardizedData) / (size(Data, 1) - 1);

% Method 3
% corrMatrix = corrcoef(Data);


% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

imagesc(corrMatrix)
colormap(jet(256))

pbaspect([1 1 1])

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Compute PCA

% Perform eigen decomposition of the covariance matrix
[eigenvectors, eigenvalues] = eig(corrMatrix);

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


%% Plot loadings (eigenvectors, principal components)
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,1,1)
stem(eigenvectors(:, 1), 'Color', [0 0 0], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 50)
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


subplot(2,1,2)
stem(eigenvectors(:, 2), 'Color', [0 0 0], 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 50)
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Sign Correction
for i = 1:size(eigenvectors, 2)
    [~, I] = max(abs(eigenvectors(:,i)));
    if eigenvectors(I,i) < 0
        eigenvectors(:,i) = -1.*eigenvectors(:,i);
    end
end


%% Project data to the principal components (eigenvectors)


PCs = standardizedData*eigenvectors(:, 1:2);

% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

subplot(5, 1, 1:3)
imagesc(Data(1:100, :)')
colormap(flip(gray(256)))
xlabel('Sample (e.g. Time)')
ylabel('Variable (e.g. Neuron/Gene)')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])


subplot(5, 1, 4)
plot(PCs(1:100, 1), 'LineWidth', 3, 'Color', [0 0 0])
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])

subplot(5, 1, 5)
plot(PCs(1:100, 2), 'LineWidth', 3, 'Color', [0 0 0])
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])

title('Activation Plot')

%% Plot trajectories
figure('units','normalized','outerposition',[0 0 1 1])

plot(PCs(:, 1), PCs(:, 2), 'LineWidth', 3, 'Color', [0 0 0])


%% Lower dimensional subspaces

PCs = standardizedData*eigenvectors(:, 1:3);

figure('units','normalized','outerposition',[0 0 1 1])

plot3(PCs(:, 1), PCs(:, 2), PCs(:, 3), 'LineWidth', 3, 'Color', [0 0 0]), hold on

% Fit a plane to the trajectory
coefficients = pca(PCs(:, 1:3));
normal_vector = coefficients(:, 3); % Normal vector to the plane
[gridX, gridY] = meshgrid(min(PCs(:,1)):0.1:max(PCs(:,1)), min(PCs(:,2)):0.1:max(PCs(:,2)));
gridZ = -(normal_vector(1) * gridX + normal_vector(2) * gridY) / normal_vector(3);

% Plot the plane
surf(gridX, gridY, gridZ, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'cyan')

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('3D Trajectory and Fitted Plane')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
grid on
hold off

%%
standardizedData = zscore(Data);
[coeff,score,latent,~,explained] = pca(standardizedData);

[coeff,score,latent,~,explained] = pca(standardizedData, 'Algorithm', 'eig');
% Algorithm: eig, svd, als

% Coeff: loadings, eigenvectors
% Score: Principal component scores are the representations of X in the principal component space.
% latent: Principal component variances, namely the eigenvalues of the covariance matrix of X,
% returned as a numeric column vector of length k.
% explained : Percentage of the total variance explained by each principal component,
% returned as a numeric column vector of length k.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, clc

%% Simulate a dataset
Data = zeros(1000, 9);
Data(randi(1000, [1 100]), 1:3) = 1;

Data(randi(1000, [1 100]), 4:6) = 1;
Data(Data(:, 1) == 1, 4:6) = 0;


Data(randi(1000, [1 100]), 7) = 1;
Data(randi(1000, [1 100]), 8) = 1;
Data(randi(1000, [1 100]), 9) = 1;


% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

imagesc(Data')
colormap(flip(gray(256)))
xlabel('Sample (e.g. Time)')
ylabel('Variable (e.g. Neuron/Gene)')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Compute correlation matrix

% Method 1
standardizedData = (Data - mean(Data))./std(Data);
corrMatrix = (standardizedData' * standardizedData) / (size(Data, 1) - 1);


% Method 2
% standardizedData = zscore(standardizedData);
% corrMatrix = cov(standardizedData) / (size(Data, 1) - 1);

% Method 3
% corrMatrix = corrcoef(Data);


% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

imagesc(corrMatrix)
colormap(jet(256))

pbaspect([1 1 1])

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Compute PCA


% Perform eigen decomposition of the covariance matrix
[eigenvectors, eigenvalues] = eig(corrMatrix);

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


%% Sign Correction
for i = 1:size(eigenvectors, 2)
    [~, I] = max(abs(eigenvectors(:,i)));
    if eigenvectors(I,i) < 0
        eigenvectors(:,i) = -1.*eigenvectors(:,i);
    end
end


%% Project data to the principal components (eigenvectors)


PCs = standardizedData*eigenvectors(:, 1:2);

% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

subplot(5, 1, 1:3)
imagesc(Data(1:100, :)')
colormap(flip(gray(256)))
xlabel('Sample (e.g. Time)')
ylabel('Variable (e.g. Neuron/Gene)')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])


subplot(5, 1, 4)
plot(PCs(1:100, 1), 'LineWidth', 3, 'Color', [0 0 0])
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])

subplot(5, 1, 5)
plot(PCs(1:100, 2), 'LineWidth', 3, 'Color', [0 0 0])
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
xlim([1 100])

title('Activation Plot')

%% Plot trajectories
figure('units','normalized','outerposition',[0 0 1 1])

plot(PCs(:, 1), PCs(:, 2), 'LineWidth', 3, 'Color', [0 0 0])


%% Lower dimensional subspaces

PCs = standardizedData*eigenvectors(:, 1:3);

figure('units','normalized','outerposition',[0 0 1 1])

plot3(PCs(:, 1), PCs(:, 2), PCs(:, 3), 'LineWidth', 3, 'Color', [0 0 0])

% Fit a plane to the trajectory
coefficients = pca(PCs(:, 1:3));
normal_vector = coefficients(:, 3); % Normal vector to the plane
[gridX, gridY] = meshgrid(min(PCs(:,1)):0.1:max(PCs(:,1)), min(PCs(:,2)):0.1:max(PCs(:,2)));
gridZ = -(normal_vector(1) * gridX + normal_vector(2) * gridY) / normal_vector(3);

% Plot the plane
surf(gridX, gridY, gridZ, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'cyan')

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('3D Trajectory and Fitted Plane')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
grid on
hold off

