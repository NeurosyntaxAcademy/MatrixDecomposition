% PART 1 - https://www.crowdcast.io/e/sara-sollas-world-wide, 15:00-25:00
% https://github.com/jbakermans/ManifoldTutorial
% Load synthetic data and find low-dimensional manifold.
% Project data on neural modes; see how latent trajectory generates data.
% Jacob Bakermans, February 2021
close all;
clear all;

%% Section 1: Load data
% In this section we'll generate some synthetic neural data for 3 neurons,
% and plot the data over time and as a trajectory in neuron activity space

% Generate synthetic data
X_D = generate_data_1(); % Data matrix: N_neurons x N_timebins
X_D = smoothdata(X_D', 'gaussian', 20)';

% Plot data: firing rate for each neuron
figure(); 
% First subplot: neuron firing rate through time
subplot(1,2,1);
h = plot(X_D', 'LineWidth', 3);
Color = [0 0 0; 0.5 0.5 0.5; 0.8 0.8 0.8];
for i = 1:3
    h(i).Color = Color(i, :);
end

% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3', 'Location', 'southwest');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Neuron activity through time');
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)
pbaspect([2.5 1 1])


% Second subplot: trajectory in neural space
subplot(1,2,2);

% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'Marker', '.', 'MarkerSize', 20,...
    'LineWidth', 3, 'Color', [0 0 0]);

% Set plot layout properties
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
view(75, 30);
grid on;
title('Trajectory in neural space');

%% Section 2: Do PCA to get neural modes/principal components 
% In this section we'll find the flat subspace/linear manifold in the data by PCA.
% This manifold is spanned by the first two principal components, or neural modes.

% A PCA is an eigendecomposition of the covariance matrix of a series of datapoints
% Remove mean firing rate for every neuron
X_D = X_D - repmat(mean(X_D, 2), [1, size(X_D,2)]);
% Calculate covariance matrix between neurons
cov_dat = X_D*X_D';
% Do eigenvalue decomposition of covariance
[V, D] = eig(cov_dat);

% V has eigenvectors in columns: these are the principal components, or neural modes
% D has corresponding eigenvalues on diagonal
% The higher the eigenvalue, the more variance of activity is explained by that neural mode
% Therefore it's useful to sort modes by eigenvalue: highest eigenvalues first

% Sort eigenvectors (columns) in V by descending eigenvalue (diagonal) in D
V = sortrows([diag(D) V'],'descend');
D = diag(V(:,1));
V = V(:,2:end)';

% Plot neural trajectory and principal components
figure();
hold on;
% Plot data
plot3(X_D(1,:), X_D(2,:), X_D(3,:), 'Marker', '.', 'MarkerSize', 20,...
    'LineWidth', 3, 'Color', [0 0 0]);

% Plot principal components
for currDir = 1:3
    quiver3(0, 0, 0, max(abs(X_D(:)))*V(1, currDir), ...
        max(abs(X_D(:)))*V(2, currDir), max(abs(X_D(:)))*V(3, currDir), 0, ...
        'LineWidth', 4);
end

% Plot plane spanned by first two principal components
fmesh(@(s,t)V(1,1)*s+V(1,2)*t, @(s,t)V(2,1)*s+V(2,2)*t, @(s,t)V(3,1)*s+V(3,2)*t, ...
    [-1, 1])
alpha(0.5);
hold off;
% Set plot layout properties
legend('Data', 'Principal component 1',...
    'Principal component 2', 'Principal component 3', 'Manifold/subpace');
xlabel('Neuron 1 activity');
ylabel('Neuron 2 activity');
zlabel('Neuron 3 activity');
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
view(75, 30);
grid on;

%% Section 3: Go from neuron activity space to neural mode space
% In this section, instead of plotting the trajectory in neuron activity space, we'll move
% the trajectory to the neural mode space: the space spanned by the neural modes
% (u and v on the slides). We'll get what Solla calls latent rates. You could think about
% latent rates as generating the measured neuron activity by linear combination.

% The first two neural modes span the manifold. Select the corresponding columns from V
V_tilde = V(:, [1 2]);

% Sign Correction
for i = 1:size(V_tilde, 2)
    [~, I] = max(abs(V_tilde(:,i)));
    if V_tilde(I,i) < 0
        V_tilde(:,i) = -1.*V_tilde(:,i);
    end
end


% Now move trajectory in neuron activity space to neural mode space to get latent rates
L = V_tilde' * X_D;

% Plot trajectory in neural mode space
figure(); 
plot(L(1,:), L(2,:), '-x');
axis square;
xlabel('Neural mode 1 (u)');
ylabel('Neural mode 2 (v)');
title('Latent variables: data projected on manifold, in neural mode space');

% Use linear combinations of the latent trajectory to recover original activity
X_gen = V_tilde * L;

% Plot neuron activity through time, and compare to latent activity and recovered activity
figure(); 
% First subplot: original data
subplot(3,1,1)
h = plot(X_D', 'LineWidth', 3);
Color = [0 0 0; 0.5 0.5 0.5; 0.8 0.8 0.8];
for i = 1:3
    h(i).Color = Color(i, :);
end
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Neuron activity through time');


% Second subplot: latent activity
subplot(3,1,2)
h = plot(L', 'LineWidth', 3);
Color = [0 0 0; 0.5 0.5 0.5; 0.8 0.8 0.8];
for i = 1:2
    h(i).Color = Color(i, :);
end
% Set plot layout properties
legend('Latent unit 1', 'Latent unit 2');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Latent activity through time');
% Third subplot: recovered original data from latent units
subplot(3,1,3)
h = plot(X_gen', 'LineWidth', 3);
Color = [0 0 0; 0.5 0.5 0.5; 0.8 0.8 0.8];
for i = 1:3
    h(i).Color = Color(i, :);
end
% Set plot layout properties
legend('Neuron 1', 'Neuron 2', 'Neuron 3');
ylabel('Normalised spike counts (1)');
xlabel('Time (bins)');
title('Recovered activity through time');


