% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy


% https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006643
% https://journals.physiology.org/doi/epdf/10.1152/jn.00254.2013
% https://github.com/BrainDynamicsUSYD/NeuroPattToolbox

clear all, clc

% load(fullfile('C:\Users\abbasps\Documents\Matrix Decomposition\NeuroPattToolbox-master\NeuroPattToolbox-master', 'sampleData_marmosetLFP.mat'))
LFP = squeeze(sampleLFP(:, :, :, 1));
Fs = sampleFs;

%%
figure

plot(squeeze(LFP(1, 5, :))', 'Color', [0 0 0], 'LineWidth', 2)
title('LFP signal from 1 example channel');
xlabel('Time');
ylabel('Amplitude');
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 20)
xlim([0 2000])


figure
for i = 1:20
    
    subplot(5, 4, i)
    imagesc(LFP(:, :, i*50))
    cmap = perpl_RedBlueColormap;
    colormap(cmap)
    caxis([-2 2]*10e-5)
    pbaspect([1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out')

end


%%
% Original size = 10 * 10 * 2035 --> reshaped size = 100 * 2035

reshapeLFP = reshape(LFP, [], size(LFP, 3))';

[U,S,V] = svd(reshapeLFP, 'econ');
Eigelvalue = sort(diag(S), 'descend'); Eigelvalue = Eigelvalue.^2;
ExpVar = cumsum(Eigelvalue/sum(Eigelvalue));

%% Explained Variance
figure('units','normalized','outerposition',[0 0 1 1])

plot(ExpVar, 'LineWidth', 3, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 100)

title('Scree Plot - Explained Variance')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)

%% Reconstruction

% Extract the dominant spatiotemporal patterns
dominantSpatialPatterns = V(:, 1:10); % First 10 spatial components
singularValues = diag(S);
dominantTemporalPatterns = U(:, 1:10); % First 10 temporal components

% Plot the dominant spatial patterns
figure;
for i = 1:10
    subplot(4, 3, i);
    imagesc(reshape(dominantSpatialPatterns(:, i), 10, 10));
    colormap(jet);
    colorbar;
    title(['Dominant Spatial Pattern ', num2str(i)]);
    axis equal tight;
end

% Plot the dominant temporal patterns
figure;
for i = 1:3
    subplot(3, 1, i);
    plot(dominantTemporalPatterns(:, i), 'LineWidth', 2);
    title(['Dominant Temporal Pattern ', num2str(i)]);
    xlabel('Time');
    ylabel('Amplitude');
    grid on;
end

%%

number_of_iterations = 100;

ICAPatterns =...
    fast_ica(reshapeLFP, 3, number_of_iterations);

prjs =  reshapeLFP * ICAPatterns;

plot(prjs)








