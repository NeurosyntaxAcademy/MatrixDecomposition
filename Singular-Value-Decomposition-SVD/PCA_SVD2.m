% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy


%% Step 1: Create the dataset
close all, clc

% Generate time vector
Fs = 1000; % Sampling frequency
duration = 10; % Duration of the signal in seconds

t = [0:1/Fs:duration]';

% Generate three main oscillations
osc1 = sin(2 * pi * 20 * t);
osc2 = cos(2 * pi * 5 * t);
osc3 = sin(2 * pi * 1 * t);

% Combine oscillations into a matrix
oscillations = [osc1, osc2, osc3];

figure,
h = plot(oscillations(1:500, 1:3), 'LineWidth', 5)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Latent Space');
xlabel('Time');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

%%


tmp = 3*osc1+0*osc2+5*osc3;

h = plot(tmp(1:500, 1), 'LineWidth', 5);




%% Project to higher dimensions using a linear function and add noise

% Random linear transformation matrix (10x3)
transformation_matrix = randn(10, 3);

% Project oscillations to 10 dimensions
X_projected = oscillations * transformation_matrix';


figure,

h = plot(X_projected(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(jet(size(X_projected, 2))); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Noisy high dimensional data');
xlabel('Time');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)



%%
% Add Gaussian noise
noise = randn(size(X_projected)) * 0.1;
X_noisy = X_projected + noise;

figure,

h = plot(X_noisy(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(jet(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Noisy high dimensional data');
xlabel('Time');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)


%% Step 2: Perform PCA using eig function

% Standardize the data matrix
X_std = zscore(X_noisy);

% Compute the covariance matrix
% n = size(X_std, 1);
% C = cov(X_std) / (n - 1);

% Perform eigen decomposition
[eigvectors, eigvalues] = eig(X_std'*X_std);

% Extract the eigenvalues from the diagonal matrix
eigvalues = diag(eigvalues);

% Sort the eigenvalues and eigenvectors in descending order
[eigvalues, idx] = sort(eigvalues, 'descend');
eigvectors = eigvectors(:, idx);

%% Step 3: Perform PCA using svd function

% Perform Singular Value Decomposition
[U, S, V] = svd(X_std, 'econ');

% Singular values (diagonal elements of S)
singValues = diag(S);

% Principal components
principal_components_svd = V;


%% Step 4: Compare the results from eig and svd

% Relationship between eigenvalues and singular values
% Eigenvalues of the covariance matrix are the squared singular values divided by (n-1)
n = size(X_std, 1);
eigvalues_from_svd = (singValues.^2); %/ (n - 1)

% Plot eigenvalues and singular values for comparison
figure;
subplot(1, 2, 1);
bar(eigvalues, 'FaceColor', [0 0 0]);
title('Eigenvalues from eig');
xlabel('Component');
ylabel('Eigenvalue');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

subplot(1, 2, 2);
bar(eigvalues_from_svd, 'FaceColor', [0 0 0]), hold on
% bar(singValues, 'FaceColor', [0.5 0.5 0.5])
title('Eigenvalues derived from SVD');
xlabel('Component');
ylabel('Eigenvalue');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)


% Add inset to the second subplot
axes('Position', [.7 .6 .25 .25]);
box on;
bar(singValues, 'FaceColor', [0.5 0.5 0.5]); pbaspect([1 1 1])
title('Singular Values');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 15);
xlabel('Component');
ylabel('Singular Value');

%%
% V = eigenvalue decomposition of A'A (e.g. X_std'*X_std)

figure;
subplot(1, 2, 1);
stem(V(:, 1), 'filled', 'Color', [0 0 0], 'LineWidth', 3);
title('V from SVD');
xlabel('Component');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

subplot(1, 2, 2);
stem(-eigvectors(:, 1), 'filled', 'Color', [0 0 0], 'LineWidth', 3);
title('Eigenvectors [A''A]');
xlabel('Component');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)

%% Reconstruction by dimension reduction

Dred = X_std * V(:, 1:3);

figure,
subplot(1, 2, 1);

h = plot(oscillations(1:500, :), 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Dimensionality reduction');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])

subplot(1, 2, 2);

h = plot(Dred(1:500, 2), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:1);
title('Dimensionality reduction');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])

%% Spectral Decomposition

ov = 0.5;   % overlap
winlen = 1024; % winlen
nff = 100;
freqBins = logspace(log10(1), log10(50), nff); %linspace(1, 150, nff); %logspace(log10(1), log10(200), nff);
fs = 1000;

% Spectral decomposition of original 3 oscillations
[pxx, f1] = pwelch(oscillations,winlen,floor(ov*winlen), freqBins, fs);


% Spectral decomposition of reduced 3-component data
[pxx_red, f_red] = pwelch(Dred,winlen,floor(ov*winlen), freqBins, fs);


% Plotting the results

figure;

% Original oscillations
subplot(2, 1, 1);

h = plot(f1, pxx, 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Spectral features - Original');
xlabel('Frequency [Hz]');
ylabel('Power');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 50])


% Original oscillations
subplot(2, 1, 2);

h = plot(f1, pxx_red, 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Spectral features - PCA');
xlabel('Frequency [Hz]');
ylabel('Power');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 50])


%% Independent Component Analysis

number_of_iterations = 100;

ICAPatterns =...
    fast_ica(X_std', 3, number_of_iterations);

prjs = X_std * ICAPatterns;


figure,
subplot(1, 2, 1);

h = plot(oscillations(1:500, :), 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Dimensionality reduction');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])

subplot(1, 2, 2);

h = plot(prjs(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('ICA Pattern');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])


%% Filtering

[b a] = butter(3, 30/(fs/2));
fltprjs = zeros(size(prjs));
for i =1:3
    fltprjs(:, i) = filtfilt(b, a, prjs(:, i));
end


figure,
subplot(1, 2, 1);

h = plot(oscillations(1:500, :), 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Dimensionality reduction');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])

subplot(1, 2, 2);

h = plot(fltprjs(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('ICA Pattern - Filtered');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])


%% Spectral Decomposition of ICA Patterns

% Spectral decomposition of reduced 3-component data
[pxx_red, f_red] = pwelch(prjs,winlen,floor(ov*winlen), freqBins, fs);


% Plotting the results

figure;

% Original oscillations
subplot(2, 1, 1);

h = plot(f1, pxx, 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Spectral features - Original');
xlabel('Frequency [Hz]');
ylabel('Power');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 50])


% Original oscillations
subplot(2, 1, 2);

h = plot(f1, pxx_red, 'LineWidth', 3)
% Set the colormap to gray
colors = colormap(turbo(3)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Spectral features - ICA');
xlabel('Frequency [Hz]');
ylabel('Power');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 50])


%%

r = 1;
LowRankData = U(:, 1:r)*S(1:r,1:r)*V(:, 1:r)';


figure,
subplot(1, 2, 1);

h = plot(X_noisy(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(turbo(10)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Dimensionality reduction');
xlabel('Noisy data');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])

subplot(1, 2, 2);

h = plot(LowRankData(1:500, :), 'LineWidth', 1)
% Set the colormap to gray
colors = colormap(turbo(10)); % Get 3 shades of gray
% Apply the colors to the plots using arrayfun
arrayfun(@(i) set(h(i), 'Color', colors(i, :)), 1:3);
title('Low Rank Approximation');
xlabel('Observations');
ylabel('Magnitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 500])




