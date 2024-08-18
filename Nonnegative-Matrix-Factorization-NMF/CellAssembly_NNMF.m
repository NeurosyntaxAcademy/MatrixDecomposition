clear all, close all % just clearing workspace
% https://github.com/tortlab/Cell-Assembly-Detection
% Detecting cell assemblies in large neuronal populations, V
% Lopes-dos-Santos 2013
% Neuronal assembly detection and cell membership specification by principal component analysis
% Lopes-dos-Santos 2011

%% Simulate a matrix of firing rates

% define mean firing rate
Network_opts.meanspikebin = 1;

% define number of neurons
Network_opts.nneurons = 100 ;

% define number of bins
Network_opts.nbins = 10000;

% above define assembly membership

% line below sets neurons 1,2,3 and 4 to be in assembly 1
Assembly_opts.assembly_neurons{1} = [1:40]; 
% line below sets neurons 5,6 and 7 to be in assembly 2
Assembly_opts.assembly_neurons{2} = [20:50]; 
% Assembly_opts.assembly_neurons{3} = [50:55]; 

% defines number of activation bins
Assembly_opts.number_of_activations = 300;

% defines mean rate in activation bins
Assembly_opts.meanspikerate_activations = 30;

% running the function
Activitymatrix = toy_simulation(Network_opts,Assembly_opts);


%% Plot spiking activity

fh = figure(1),clf
imagesc(Activitymatrix(:, 1:200))
colormap(flip(gray(256)));
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])
pbaspect([2 1 1])
xlabel('Time')
ylabel('Neuron#')


%% Plot correlation matrix

correlationmat = corr(Activitymatrix');
fh = figure(2),clf
imagesc(correlationmat)
colormap(perpl_RedBlueColormap); colorbar()
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
caxis([-1 1])
pbaspect([1 1 1])
xlabel('Neuron#')
ylabel('Neuron#')

%%
Options = [];
Options = statset('MaxIter', 1000, 'TolFun', 1e-4, 'UseParallel', false, 'TolX', 1e-4);

numComponents = 2; % Set the number of components for the factorization
[W, H] = nnmf(Activitymatrix', numComponents, 'options', Options);

%%
Time = 1:200;

% Plot the original image
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(5, 5, [1:3 6:8 11:13]);
imagesc(Activitymatrix(:, Time));
colormap(flip(gray));
title('Spiking Matrix');
ylabel('Neuron')
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'xtick', [], 'ytick', [1 50 100])

% Plot the H matrix
subplot(5, 5, [4:5 9:10 14:15]);
imagesc(H');
colormap(flip(gray));
colorbar;
title('W Matrix');
xlabel('Components');
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'xtick', 1:2, 'ytick', [1 50 100])
grid on

% Plot the W matrix
subplot(5, 5, [16:18 21:23]);
imagesc(W(Time, :)');
colormap(flip(gray));
title('H Matrix');
xlabel('Time');
ylabel('Component');
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'ytick', 1:2)
grid on

%% Noise reduction

% ProjAct = W*H;
% 
% fh = figure(1),clf
% imagesc(ProjAct(:, Time))
% colormap(flip(gray(256)));
% set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% % caxis([0 1])
% pbaspect([2 1 1])
% xlabel('Time')
% ylabel('Neuron#')
% 







