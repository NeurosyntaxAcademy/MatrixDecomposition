% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy

%%

clear all, clc

%% Simulate a dataset
% Data = randi(100, [1000, 1]);
% Data(:, 2) = 2*Data(:, 1); %+rand([1000, 1]) * 5


%% Dataset 2
% Parameters
Correlation = 0.8;

% Generate two independent random variables
x1 = randn(1000, 1);
x2 = randn(1000, 1);

% Combine them with the desired correlation
x2 = Correlation * x1 + sqrt(1 - Correlation^2) * x2;

% Create the data matrix
Data = [x1, x2];



%% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

scatter(Data(:, 1), Data(:, 2), 20, [0 0 0], 'Filled' )

xlabel('Va1')
ylabel('Va2')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)



%% Compute correlation matrix

% Method 1
% standardizedData = (Data - mean(Data))./std(Data);
% corrMatrix = (standardizedData' * standardizedData) / (size(Data, 1) - 1);

% Method 2
% standardizedData = zscore(standardizedData);
% corrMatrix = cov(standardizedData) / (size(Data, 1) - 1);

% Method 3
standardizedData = zscore(Data);
corrMatrix = corrcoef(standardizedData);


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

plot(eigenvalues/sum(eigenvalues)*100, 'LineWidth', 3, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 100)

title('Scree Plot - Explained Variance')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
ylabel('EV (%)')
xlabel('Principal Component')


%%


figure('units','normalized','outerposition',[0 0 1 1])

scatter(standardizedData(:, 1), standardizedData(:, 2), 20, [0.5 0.5 0.5], 'Filled' ), hold on

for i = 1:2
    hold on
    quiver(0, 0, eigenvectors(1, i), eigenvectors(2, i), ...
        'LineWidth', 2, 'MaxHeadSize', 2, 'AutoScale', 'on', 'AutoScaleFactor', eigenvalues(i));
end


xlabel('Va1')
ylabel('Va2')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)


%% Project data to the principal components (eigenvectors)

PCs = standardizedData*eigenvectors(:, 1:2);

% Plot data
figure('units','normalized','outerposition',[0 0 1 1])

scatter(PCs(:, 1), PCs(:, 2), 20, [0 0 0], 'Filled' )

xlabel('Va1')
ylabel('Va2')

set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)

xlim([min(PCs(:, 1))-1 max(PCs(:, 1))+1])
ylim([min(PCs(:, 1))-1 max(PCs(:, 1))+1])







