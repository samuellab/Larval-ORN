% demonstrate how to pick one model from the other
% model A, each response curve's three parameters are all different
% model B, they share the same slop
% figure out the best way to estimate the parameters for each curve
% compare the R^2 of the fitting

% run after the scirpt 'FitSingleORNData.m', which saved the data and some
% fitting resutls

%%
close all; clear; clc;
warning('off','all');
diary(fullfile('.', 'AnalysisResults', 'modelcomparison.txt')); 
diary on;

%% load data
fName = fullfile('AnalysisResults', 'fitSingleORNdataResults.mat');
load(fName);

% format data
dataX = [];	dataY = [];  coeffA = [];
for i = 1:2
    for j= 1:2
        dataY = [dataY; input.dff{i,j}];
        [trailNum, ~] = size(input.dff{i,j});
        dataX = [dataX; log10(repmat(input.concList{i,j}, trailNum, 1))];
        
        coeffA = [coeffA; results.fitCoeffIdv{i,j}];
    end
end

% model funciton
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

%% model A
% the best estimated parameters are the ones best for each individual
% curve, which was fitted before. So, here, calculate the R^2.
rSquareA = EnsmbleRSquare(dataX, dataY, coeffA, hillEq);
title(['Model A, R^2=', num2str(rSquareA)]); 
xlabel('$$y$$','Interpreter','Latex'); ylabel('$$\hat{y}$$','Interpreter','Latex');

disp('----------MODEL A----------')
disp(['R^2 = ', num2str(rSquareA)]);

%% model B
slop0 = median(coeffA(:, 2)); ampVec0 = coeffA(:, 1);	kdVec0  = coeffA(:, 3);
[slop, ampVec, kdVec, ~] = EnsembleMiniSearch(dataX, dataY, hillEq, slop0, ampVec0, kdVec0);

coeffB = [ampVec repmat(slop, length(ampVec), 1), kdVec];
rSquareB = EnsmbleRSquare(dataX, dataY, coeffB, hillEq); 
title(['Model B, R^2=', num2str(rSquareB)]); 
xlabel('$$y$$','Interpreter','Latex'); ylabel('$$\hat{y}$$','Interpreter','Latex');

disp('----------MODEL B----------')
disp(['R^2 = ', num2str(rSquareB)]);
disp(['Slop = ', num2str(slop)]);

%% compare the parameters
disp('----------COMPARE PARAMETERS----------');
fprintf('%5s\t%-5s\t%-5s\t%-5s\t%-5s\t%-5s\t%5s\t\n', '#', 'Amp_A', 'Amp_B', 'Slop_A', 'Slop_B','EC50_A', 'EC50_B');
for i = 1:length(ampVec)
    fprintf('%5s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', num2str(i), ...
        coeffA(i, 1), coeffB(i, 1), coeffA(i, 2), coeffB(i, 2), ...
        coeffA(i, 3), coeffB(i, 3));
end

%% davedata
mcResults.dataX = dataX;    mcResults.dataY = dataY;
mcResults.coeffA = coeffA;  mcResults.coeffB = coeffB;
mcResults.r2A = rSquareA;   mcResults.r2B = rSquareB;

save(fullfile('.', 'AnalysisResults', 'modelComparisonResults.mat'), 'mcResults');
diary off;

% save all figures
FolderName = fullfile('.', 'AnalysisResults');   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = ['modelComparison', num2str(get(FigHandle, 'Number')), '.fig'];
  savefig(FigHandle, fullfile(FolderName, FigName));
end