%% load the fitted EC50 matrix
clear; clc;
warning('off','all');
% diary(fullfile('.', 'AnalysisResults', 'AnalyzeEC50Matrix_log.txt')); 
% diary on;

load(fullfile('.', 'AnalysisResults', 'FitDoseResponse.mat'));

% add two super senstive pairs, remove after measurement 
cMatrix(23, 12) = -8.8;   %2-heptanone, 85c
cMatrix(24, 16) = -9.2;  %methyl salicylate, 22c 

% load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
load(fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat'));

% %%
% figure;
% imagesc(-cMatrix);

%% distribution of the value of 1/EC50
% pool the non-NaN elements
% ec50Pool = cMatrix(~isnan(cMatrix));
ec50Pool = nonzeros(cMatrix);

sData = 1./(10.^ec50Pool);

% log-log plot of the CDF of the values of 1/EC50
% CDF, count from the large to small
sDataS = sort(sData);
yPlot = (length(sDataS):-1:1)/length(cMatrix(:));
figure; loglog(sDataS ,yPlot, 'o' ); 
xlabel('1/EC50'); ylabel('CDF'); title('CDF(1/EC50)');

%% fit and verify the power-law distribution
% Code from the paper of 'Power-law Distributions in Empirical Data',
% downloaded from http://tuvalu.santafe.edu/~aaronc/powerlaws/

% add the downloadeded code to current work directory
addpath(fullfile('.', '3rdPartyCodes', 'powerlaws_full_v0.0.10-2012-01-17'));

% fit the power law function
[alpha, xmin, L] = plfit( sData );
% visualize the fiting results
plplot(sData, xmin, alpha);
% estimates the uncertainty in the estimated power-law parameters.
[p,gof] = plpva(sData, xmin, 'silent');

% print the fitting resutls
disp('----------Fit Distribution Function of 1/EC50 Data----------');
disp('Reference: http://tuvalu.santafe.medu/~aaronc/powerlaws/');
disp('PDF: p(x) = x^-alpha, for x >= xin');
fprintf('%25s: alpha = %.2f\n', 'Scaling exponent', alpha);
fprintf('%25s: xmin = %.2e\n', 'Lower bound', xmin);
fprintf('%25s: L = %.2e\n','Log-likelihood', L);
fprintf('%25s: p = %.2f\n', 'p-Value', p);
fprintf('(Quote from the paper) If the resulting p-value is greater than 0.1,\n the power law is a plausible hypothesis for the data, otherwise it is rejected.\n')

%% apply PCA on the EC50 matrix
% replace NaN in the EC50 matrix with 0
% cMatrix(isnan(cMatrix)) = 0;

% PCA, using -cMatrix does not change the dimentionality, but it could make
% the load and projection match each other
[coeff1, score1, latent1,tsquared1,explained1, ~] = pca(-cMatrix);

% show the percentage of variance each PC explained
figure(); pareto(explained1);  

% plot the load/projection of/on the 1st PC.
load1st = coeff1(:, 1);
proj1st = score1(:, 1);

[load_sorted, load_index] = sort(load1st);
[proj_sorted, proj_index] = sort(proj1st);

figure; barh(load_sorted); 
set(gca,'YTick',1:length(load_index));
set(gca,'YTickLabel',ORNList(load_index));
xlabel('Load on 1st PC');

figure; barh(proj_sorted); 
set(gca,'YTick',1:length(proj_index));
set(gca,'YTickLabel',odorList(proj_index));
xlabel('Projection on 1st PC');

%% Compare odors' functional with structure
%load structure data
%load the .xlsx file of the odor strucure descriptors
filename = fullfile('.', 'data', 'Odor_Structure_Descriptors_EDragon.xlsx');
sheet = 1;
xlRange = 'A3:BLF37';
[~,~,raw] = xlsread(filename,sheet,xlRange);

infoMetric = raw(1, 5:end); %define the information list of the descriptors
m = cell2mat(raw(2:end, 5:end)); %extract the data matrix

% Select the optimized descriptor from 'A metric for odorant comparison'
optList  = [22 48 92 96 359 404	433	513	582	592	678	688	690	698	946 ...
    948	959	963	1012 1069 1110 1191 1286 1321 1331 1348 1373 1430 1528 ...
    1541 1558 1576];

inforMOpt = infoMetric(optList);
mOpt = m(:, optList);

% calculate the Pearson correlation coefficient between functional
% projection and each structure metric
rho = zeros(1, length(mOpt));
for i = 1:length(optList)
    temp = corrcoef(proj1st, mOpt(:,i));
    rho(i) = temp(1, 2);
end
%find the largest correlation coefficient
[rhoMax, rhoMaxInd] = max(rho);

% calculate the linear fitting
xdata = proj1st;
ydata = mOpt(:, rhoMaxInd);
pPoly = polyfit(xdata, ydata, 1); % Linear fit of xdata vs ydata
linePointsX = [min(xdata) max(xdata)]; % find left and right x values
linePointsY = polyval(pPoly, [min(xdata),max(xdata)]); % find y values

% correlation plot
figure; 
plot(xdata, ydata, 'ko'); hold on
plot(linePointsX,linePointsY,'-r');
xlabel('Projection on 1st PC');
ylabel(inforMOpt{rhoMaxInd});
title(['Correlation: ', num2str(rhoMax, 1)]);

% display the results 
disp('----------Relate Odorant Functional Data with Structural Data----------');
disp('Reference: http://www.nature.com/nmeth/journal/v5/n5/full/nmeth.1197.html');
disp('');
fprintf('Odor structure descriptor "%5s" hits the highest coefficient score with the 1st PC of odor functional data. \n', inforMOpt{rhoMaxInd});
fprintf('The Pearson correlation coefficient is %.2f. \n', rhoMax);

%% shuffle the matrix and apply PCA.
[rows, cols] = size(cMatrix);
cMShuffled = zeros(rows, cols);

N = 1000;
explainedMatrix = zeros(cols , N);

for k = 1:N
    for i = 1:rows
        %random shuffle the elements in each row
        temp = cMatrix(i, :);
        cMShuffled(i, :) = temp(randperm(cols));
    end
    %apply PCA
    [~, ~, ~, ~, explained, ~] = pca(-cMShuffled);
    %keep the varaince explained
    explainedMatrix(:, k) = explained;
end

%calculate the mean and std of the percentage explained
explainedAve = mean(explainedMatrix, 2);
explainedStd = std(explainedMatrix, 0, 2);

% figure;   pareto(explainedAve);	title('Shuffled Matrix');

figure;   plot(1:length(explained1), explained1, '-ok'); hold on;
errorbar(1:length(explainedAve), explainedAve, explainedStd, 's');
patch([1:length(explainedAve) fliplr(1:length(explainedAve))], ...
    [explainedAve'+explainedStd' fliplr(explainedAve'-explainedStd')],[0.7 0.7 0.7]);
xlabel('PCs'); ylabel('% of variance explained'); legend('KdMatrix', 'Shuffled');
hold off;

figure; hist(explainedMatrix(1,:)); hold on; plot(explained1(1)*ones(1, 31), 0:10:300, 'r'); 
xlabel('% of variance explained by 1st PC'); ylabel('counts');
mu1 = mean(explainedMatrix(1,:)); sigma1 = std(explainedMatrix(1,:)); 
significance = (explained1(1) - mu1)/sigma1 ; 
title([num2str(round(significance)), ' sigma']);

figure; hist(explainedMatrix(2,:)); hold on; plot(explained1(2)*ones(1, 31), 0:10:300, 'r'); 
xlabel('% of variance explained by 2nd PC'); ylabel('counts');
mu2 = mean(explainedMatrix(2,:)); sigma2 = std(explainedMatrix(1,:)); 
significance2 = (explained1(2) - mu2)/sigma2 ; 
title([num2str(round(significance2)), ' sigma']);

disp('----------PCA on the log(1/EC50) Data----------');
fprintf('1st PC accounts %.2f percentage of data variance. \n', explained1(1));
fprintf('%.0f sigma significance compare with shuffled data. \n', round(significance));

%% Save analyzed data 
save(fullfile('.', 'AnalysisResults', 'AnalyzeEC50Matrix.mat'), 'alpha', 'xmin', 'p');
% diary off;