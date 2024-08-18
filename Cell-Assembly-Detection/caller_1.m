clear all % just clearing workspace
% https://github.com/tortlab/Cell-Assembly-Detection
% Detecting cell assemblies in large neuronal populations, V
% Lopes-dos-Santos 2013
% Neuronal assembly detection and cell membership specification by principal component analysis
% Lopes-dos-Santos 2011
% Edited by Saman Abbaspoor 2024

% Refs: Detecting cell assemblies in large neuronal populations, 2013
% Neuronal assembly detection and cell membership specification by
% principal component analysis, 2011


%% Simulate a matrix of firing rates

% define mean firing rate
Network_opts.meanspikebin = 1;

% define number of neurons
Network_opts.nneurons = 100 ;

% define number of bins
Network_opts.nbins = 10000;

% above define assembly membership

% line below sets neurons 1,2,3 and 4 to be in assembly 1
Assembly_opts.assembly_neurons{1} = [1:30]; 
% line below sets neurons 5,6 and 7 to be in assembly 2
Assembly_opts.assembly_neurons{2} = [20:40]; 
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


%% Compute PCA and find patterns

opts.threshold.permutations_percentile = 95;
opts.threshold.number_of_permutations = 20;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'PCA';
AssemblyTemplates = assembly_patterns(Activitymatrix,opts);
AssemblyTemplatesPCA = AssemblyTemplates;

%% Plot loadings and interpret them

figure(3),clf
subplot(211)
stem(AssemblyTemplates(:,1), 'Color', [0 0 0], 'LineWidth', 2)
% ylim([-0.2 0.6])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'ytick', [0 0.5])


subplot(212)
stem(AssemblyTemplates(:,2), 'Color', [0 0 0], 'LineWidth', 2)
% ylim([-0.2 0.6])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'ytick', [0 0.5])
xlabel('Neuron#')
ylabel('PCA Weight')


%% Activation Vector
% Trick from: https://github.com/DavidTingley/RnR_methods/tree/master/reactivation

scoreTar    = zscore(Activitymatrix')*AssemblyTemplates;

% This trick is used to get rid of the diagonal term in react. strength
tmp     = (zscore(Activitymatrix').^2)*(AssemblyTemplates.^2);
ActivationStrength       = scoreTar.^2 - tmp;

ActivationStrength = ActivationStrength';


%% Plot spiking activity and cell assembly activation matrix

fh = figure(1),clf
subplot(5,1,[1 2 3])
Time = 1:200;
imagesc(Activitymatrix(:, Time)), xlim([min(Time) max(Time)])
colormap(flip(gray(256)));
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])
ylabel('Neuron#')


subplot(5,1,[4])
Time = 1:200;
plot(ActivationStrength(1, Time), 'Color', [137 98 121]/255, 'LineWidth', 3), xlim([min(Time) max(Time)])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])


subplot(5,1,[5])
Time = 1:200;
plot(ActivationStrength(2, Time), 'Color', [173 178 211]/255, 'LineWidth', 3), xlim([min(Time) max(Time)])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])
xlabel('Time')
ylabel('Activation Strength')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform ICA
opts.threshold.permutations_percentile = 95;
opts.threshold.number_of_permutations  = 20;
opts.Patterns.number_of_iterations     = 1000;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
AssemblyTemplates = assembly_patterns(Activitymatrix(:, 5001:end),opts);
AssemblyTemplatesICA = AssemblyTemplates;


%% Sign Correction
% PCA/ICA methods are exact up to the sign of the loadings.
% Convention: the arbitrary signs of the weights were set so that the highest absolute weight was positive.
% Refs: Flexible communication between cell assemblies and 'reader'neurons,
% BioRxiv, 2022

for i = 1:size(AssemblyTemplates, 2)
    [~, I] = max(abs(AssemblyTemplates(:,i)));
    if AssemblyTemplates(I,i) < 0
        AssemblyTemplates(:,i) = -1.*AssemblyTemplates(:,i);
    end
end

%% Plot loadings and interpret them

figure(4),clf
subplot(211)
stem(AssemblyTemplates(:,1), 'Color', [0 0 0], 'LineWidth', 2)
% ylim([-0.2 0.6])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'ytick', [0 0.5])


subplot(212)
stem(AssemblyTemplates(:,2), 'Color', [0 0 0], 'LineWidth', 2)
% ylim([-0.2 0.6])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20, 'ytick', [0 0.5])
xlabel('Neuron#')
ylabel('PCA Weight')


%% Activation Vector

scoreTar    = zscore(Activitymatrix')*AssemblyTemplates;

% This trick is used to get rid of the diagonal term in react. strength
tmp     = (zscore(Activitymatrix').^2)*(AssemblyTemplates.^2);
ActivationStrength       = scoreTar.^2 - tmp;

ActivationStrength = ActivationStrength';


%% Plot spiking activity and cell assembly activation matrix

fh = figure(1),clf
subplot(5,1,[1 2 3])
Time = 1:200;
imagesc(Activitymatrix(:, Time)), xlim([min(Time) max(Time)])
colormap(flip(gray(256)));
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])
ylabel('Neuron#')


subplot(5,1,[4])
Time = 1:200;
plot(ActivationStrength(1, Time), 'Color', [137 98 121]/255, 'LineWidth', 3), xlim([min(Time) max(Time)])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])


subplot(5,1,[5])
Time = 1:200;
plot(ActivationStrength(2, Time), 'Color', [173 178 211]/255, 'LineWidth', 3), xlim([min(Time) max(Time)])
set(gca, 'Box', 'off', 'tickdir', 'out','fontsize', 20)
% caxis([0 1])
xlabel('Time')
ylabel('Activation Strength')


    