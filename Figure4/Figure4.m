clear;  clc; warning('off','all');

diary(fullfile('.', 'results', 'Figure4_log.txt')); 
diary on;

%% load the fitted log_10(EC50) matrix
load(fullfile('.', 'data', 'log10EC50.mat'));
cMatrix = log10EC50;

%% apply PCA on the EC50 matrix
% replace NaN in the EC50 matrix with 0
cMatrix(isnan(cMatrix)) = 0;

% PCA. We perform PCA on the -cMatrix -> this does not change the 
% dimentionality, ut allows coeff(loads) and score(projection) axes to
% match such that  both axes go from aromatic to aliphatic
[coeff1, score1, latent1,tsquared1,explained1, ~] = pca(-cMatrix);

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

%% Compare function with structure of odorants
%load structure data
%load the .xlsx file of the odor structure descriptors
filename = fullfile('.', 'data', 'Odor_Structure_Descriptors_EDragon.xlsx');
sheet = 1;
xlRange = 'A3:BLF37';
[~,~,raw] = xlsread(filename,sheet,xlRange);

infoMetric = raw(1, 5:end); %define the information list of the descriptors
m = cell2mat(raw(2:end, 5:end)); %extract the data matrix

% Select the optimized descriptor from paper: 'A metric for odorant comparison'
optList  = [22 48 92 96 359 404	433	513	582	592	678	688	690	698	946 ...
    948	959	963	1012 1069 1110 1191 1286 1321 1331 1348 1373 1430 1528 ...
    1541 1558 1576];

inforMOpt = infoMetric(optList);
mOpt = m(:, optList);

% calculate the Pearson correlation coefficient between functional
% projection and each structure metric
rho = zeros(1, length(mOpt));
for i = 1:length(optList)
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

% show the results 
disp('----------Relate Odorant Functional Data with Structural Data----------');
disp('Reference: http://www.nature.com/nmeth/journal/v5/n5/full/nmeth.1197.html');
disp('');
fprintf('Odor structure descriptor "%5s" hits the highest coefficient score with the 1st PC of odor functional data. \n', inforMOpt{rhoMaxInd});
fprintf('The Pearson correlation coefficient is %.1f. \n', rhoMax);

%%
diary off;