% Saman Abbaspoor 2024
% https://github.com/NeurosyntaxAcademy


%% Load image

A = imread('street1.jpg');
A = rgb2gray(A);
imshow(A)
title(['Original (',sprintf('Rank %d)',rank(double(A)))])


%% Compute SVD

[U,S,V] = svd(double(A), 'econ');
SinValue = sort(diag(S), 'descend'); Eigelvalue = SinValue.^2;
ExpVar = cumsum(Eigelvalue/sum(Eigelvalue));

%% Explained Variance

figure('units','normalized','outerposition',[0 0 1 1])

plot(ExpVar, 'LineWidth', 3, 'Color', [0 0 0], 'Marker', '.', 'MarkerSize', 50)

title('Scree Plot - Explained Variance')
set(gca, 'Box', 'off', 'tickdir', 'out', 'FontSize', 20)

%%
r = 360;
LowRimage = U(:, 1:r)*S(1:r,1:r)*V(:, 1:r)';

figure
subplot(121)
imshow(A)
title(['Original (',sprintf('Rank %d)',rank(double(A)))])


subplot(122)
imshow(uint8(LowRimage))
title(sprintf('Rank %d approximation',r))


%%
Data = A;

beta = size(Data,1) / size(Data,2);
thresh = optimal_SVHT_coef(beta, 0);
[U, S, V] = svd(double(A));
% S( S < thresh ) = 0;
% Xhat = U * S * V';

figure
SinVal = log(diag(S));
plot(1:length(SinVal), SinVal, 'Marker', '.', 'MarkerSize', 20, 'Color', [0 0 0])
hold on
TruncatedSinVal = (SinVal < thresh) .* SinVal; TruncatedSinVal(TruncatedSinVal == 0) = NaN;
plot(TruncatedSinVal, 'Marker', '.', 'MarkerSize', 20, 'Color', [1 0 0])


figure
subplot(121)
imshow(A)

subplot(122)
imshow(uint8(Xhat))
r = sum(logical(S), 'all');




