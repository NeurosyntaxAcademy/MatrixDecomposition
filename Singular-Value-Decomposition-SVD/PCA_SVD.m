% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy

%%

% Generate a random dataset (e.g., 100 samples, 5 features)
data = randn(100, 5);

% Standardize the data (zero mean, unit variance)
data = (data - mean(data)) ./ std(data);

% Perform SVD on the standardized data
[U, S, V] = svd(data); %'econ'
singularValues = diag(S);

% Perform PCA using eig on the correlation matrix
[eigenvectors, eigenvaluesMatrix] = eig(data'*data);
eigenvalues = diag(eigenvaluesMatrix);
[eigenvalues, idx] = sort(eigenvalues, 'descend'); % Sort eigenvalues in descending order
eigenvectors = eigenvectors(:, idx);

% The singular values should be equal to the square root of the eigenvalues
sqrtEigenvalues = sqrt(eigenvalues);

% Display the results
disp('Singular Values from SVD:');
disp(singularValues);

disp('Square Root of Eigenvalues from PCA:');
disp(sqrtEigenvalues);

% Plot the comparison
figure;
subplot(1, 2, 1);
plot(singularValues, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
title('Singular Values from SVD');
xlabel('Index');
ylabel('Value');
grid on;

subplot(1, 2, 2);
plot(sqrtEigenvalues, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
title('Square Root of Eigenvalues from PCA');
xlabel('Index');
ylabel('Value');
grid on;

% Superimpose both plots to visualize the equality
figure;
plot(singularValues, 'bo-', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(sqrtEigenvalues, 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
legend('Singular Values from SVD', 'Square Root of Eigenvalues from PCA');
title('Comparison of Singular Values and Square Root of Eigenvalues');
xlabel('Index');
ylabel('Value');
grid on;
hold off;

%%
figure;

imagesc(eigenvectors)

size(eigenvectors)
size(V)


