%% load the fitted EC50 matrix
clear; 
% clc;
warning('off','all');
% diary(fullfile('.', 'AnalysisResults', 'AnalyzeEC50Matrix_log.txt')); 
% diary on;

% load(fullfile('.', 'AnalysisResults', 'FitDoseResponse.mat'));

% load(fullfile('.', 'AnalysisResults', 'fitResults.mat'));

% load('\\LABNAS100\Guangwei\code\fitting\FinalFittingResult\MLE_fit_all_20180718.mat')
% cMatrix = c0MLE;

load(fullfile('.', 'AnalysisResults', 'MLEFit_Final_20180724.mat'));
cMatrix = cMatrixMLE;


% load(fullfile('.', 'data', 'AveRawDataMatrix.mat'));
load(fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat'));

% %%
% figure;
% imagesc(-cMatrix);
%% visualize the EC50 matrix

ORNOrder = [19,21,3,6,8,14,10,16,9,7,11,4,12,13,1,2,20,5,17,15,18];
odorOrder = [19,33,12,32,27,15,14,7,8,9,6,22,24,31,30,29,1,5,10,3,23,4,20,28,26,17,13,25,16,2,21,11,34,18];

[rowNum, colNum] = size(cMatrix);
newMStep1 = -cMatrix;
for i = 1:rowNum
    newMStep1(i,:) = -cMatrix(odorOrder(i),:);
end
newMStep2 = newMStep1;
for i = 1:colNum
    newMStep2(:,i) = newMStep1(:,ORNOrder(i));
end

% draw the heat map of the EC50 matrix
ec50Map = newMStep2;
figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 610, 420]);
imagesc(ec50Map); 
set(gcf, 'Position', [100 250 560 700])
set(gca, 'CLim', [0 max(ec50Map(:))]);
set(gca,'XTick',1:colNum);
set(gca,'XTickLabel', ORNList(ORNOrder));
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:rowNum);
set(gca,'YTickLabel', odorList(odorOrder));
set(gca, 'XTickLabelRotation', 45);
cmp = colormap(jet); cmp(1,:) = [0 0 0];
colormap(cmp); c = colorbar; 
c.TickLabels{1} = 'NaN'; c.Label.String = '-log10(EC50)';
title('EC50'); 


%% distribution of the value of 1/EC50
% pool the non-NaN elements
ec50Pool = cMatrix(~isnan(cMatrix));
% ec50Pool = nonzeros(cMatrix);
% ec50Pool = nonzeros(cMatrixMol);

% ec50Pool = ec50Pool(ec50Pool<-3.5);

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
[alpha, xmin, L] = plfit( sData, 'limit', 1e+5);

% visualize the fiting results
plplot(sData, xmin, alpha);
% estimates the uncertainty in the estimated power-law parameters.
[p,gof] = plpva(sData, xmin, 'silent');

% print the fitting resutls
disp('----------Fit Distribution Function of 1/EC50 Data----------');
disp('Reference: http://tuvalu.santafe.medu/~aaronc/powerlaws/');
disp('PDF: p(x) = x^-alpha, for x >= xim');
fprintf('%25s: alpha = %.2f\n', 'Scaling exponent', alpha);
fprintf('%25s: xmin = %.2e\n', 'Lower bound', xmin);
fprintf('%25s: L = %.2e\n','Log-likelihood', L);
fprintf('%25s: p = %.2f\n', 'p-Value', p);
fprintf('(Quote from the paper) If the resulting p-value is greater than 0.1,\n the power law is a plausible hypothesis for the data, otherwise it is rejected.\n')


%% randomly assemble a subset of odors calculate the distribution 
N = 1000;
nOdor = 18;

fitSubCoef = zeros(N, 3);

% disp('----------Fit 1/EC50 Distribution Using Subset of Data----------');
% fprintf('%-3s\t%-4s\t%-5s\t%-3s\t\n', '#', 'alpha', 'xmin', 'p');

parfor i = 1:N
    list = randperm(34);
    list = list(1:nOdor);
    cMatSubset = cMatrix(list, :);

    cPoolSub = cMatSubset(~isnan(cMatSubset));
    sDataSub = 1./(10.^cPoolSub);
    sDataSub = sort(sDataSub);
    
    [alphaSub, xminSub, LSub] = plfit( sDataSub, 'limit', 1e+5);
    [pSub,gofSub] = plpva(sDataSub, xminSub, 'silent');
    
    fitSubCoef(i, :) = [alphaSub, xminSub, pSub];
    
% 	fprintf('%.0f\t%.2f\t%.0f\t%.2f\t\n',i, alphaSub, xminSub, pSub);

end
%%
disp(['Count of P-vlaue greater than 0.1: ', num2str(length(find(fitSubCoef(:, 3) > 0.1)))]);
figure; histogram(fitSubCoef(:,end), 10);
idxSub = find(fitSubCoef(:, 3) > 0.1);
figure; histogram(fitSubCoef(idxSub, 1), 10); title([num2str(mean(fitSubCoef(idxSub, 1))), ' +- ', num2str(std(fitSubCoef(idxSub, 1)))])

%% apply PCA on the EC50 matrix
% replace NaN in the EC50 matrix with 0
cMatrix(isnan(cMatrix)) = 0;

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

%%
figure; stem(proj_sorted); 

figure; stem(load_sorted); 

% set(gca,'YTick',1:length(proj_index));
% set(gca,'YTickLabel',odorList(proj_index));
% xlabel('Projection on 1st PC');

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
%     temp = corrcoef(proj1st, mOpt(:,i));
%     rho(i) = temp(1, 2);
    rho(i) = corr(proj1st, mOpt(:,i));
end
%find the largest correlation coefficient
[rhoMax, rhoMaxInd] = max(abs(rho));

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
explainedStd = std(explainedMatrix, 0, 2)*2;

% figure;   pareto(explainedAve);	title('Shuffled Matrix');

figure;   plot(1:length(explained1), explained1, '-ok'); hold on;
errorbar(1:length(explainedAve), explainedAve, explainedStd, 's');
patch([1:length(explainedAve) fliplr(1:length(explainedAve))], ...
    [explainedAve'+explainedStd' fliplr(explainedAve'-explainedStd')],[0.7 0.7 0.7]);
xlabel('PCs'); ylabel('% of variance explained'); legend('KdMatrix', 'Shuffled');
hold off;
axis([0 21 ylim]);

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
% save(fullfile('.', 'AnalysisResults', 'AnalyzeEC50Matrix.mat'), 'alpha', 'xmin', 'p');
% diary off;

%% sort the -log10(EC50) for each odor
negC = -cMatrix;
negCMax = max(negC, [], 2);
[~, rowIdx] = sort(negCMax, 'descend');

figure; 
for i = 1:rowNum
    idx = rowIdx(i);
    kdVec = negC(idx, :);
    kdVec = kdVec(find(kdVec));
    [kdVec, ~] = sort(kdVec);   
    
    iShow = i;
    
    yPlot = repmat([iShow-0.3; iShow+0.3], [1, length(kdVec)]);
    xPlot = repmat(kdVec, [2, 1]);
    plot(xPlot, yPlot, 'k'); 
    hold on;
    plot([1, 10], [iShow, iShow], 'k');
%     plot([hc(idx); hc(idx)], [iShow-0.3; iShow+0.3],'r');
end
xlabel('-log_{10}EC_{50}');axis tight; 
% set(gca,'box','off','ycolor','w')
set(gca,'box','off');
set(gcf, 'Position', [200 10 560 988]);
yticks(1 : rowNum);
yticklabels(odorList(rowIdx));

