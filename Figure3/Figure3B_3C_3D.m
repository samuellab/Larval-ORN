%% load the fitted EC50 matrix
load(fullfile('.', 'results', 'MLEFit.mat'));
cMatrix = cMatrixMLE;

load(fullfile('.', 'data', 'doseResponseData.mat'));

%% Plot figure 3B
% order for visualize
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
title('log10(EC50)'); 

saveas(gcf, fullfile('results', 'figures', 'Figure3B.fig'));
%% Plot figure 3C
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
    plot([2, 10], [iShow, iShow], 'k');

end
xlabel('-log_{10}EC_{50}');axis tight; 
set(gca,'box','off');
set(gcf, 'Position', [200 10 560 988]);
yticks(1 : rowNum);
yticklabels(odorList(rowIdx));

saveas(gcf, fullfile('results', 'figures', 'Figure3C.fig'));

%% fit and verify the power-law distribution, plot figure 3D
% Code from the paper of 'Power-law Distributions in Empirical Data',
% downloaded from http://tuvalu.santafe.edu/~aaronc/powerlaws/

addpath(genpath(pwd));

% pool the non-NaN elements
ec50Pool = cMatrix(~isnan(cMatrix));
sData = 1./(10.^ec50Pool);

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

saveas(gcf, fullfile('results', 'figures', 'Figure3D.fig'));