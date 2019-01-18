% Warnning: run this code will take a very long time...
% final optimization result is saved as ./results/MLEFit.mat

%% Load the raw dose response data and initial guesses from least squares fitting
addpath(genpath(pwd));

fileName =  fullfile('data', 'Data S1.csv');

dataT = readtable(fileName);    % load Excel file

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');	%ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% mark the type of each curve
maskPlot = SetCurveType(34, 21); % manually define the type of the dose-response curves

%% pool all the saturated data
t123idx = maskPlot >= 1 & maskPlot<4;

load(fullfile('results', 'initialParameters.mat'));
maskPlot(t123idx & isnan(cMatrix)) = 0;
t1idx = maskPlot ==1;
t123idx = maskPlot >= 1 & maskPlot<4;
[sRow, sCol] = find(t123idx);

pairNum = length(sRow);

concT123 = NaN(5, 10, pairNum);
rspT123  = NaN(5, 10, pairNum);

for i = 1 : pairNum
    odIdx  = find(input.indOdor == sRow(i)); % row index of each odor
    ornIdx = sCol(i);
    
    concPool = input.concList(odIdx);       % concentration list
    dffPool  = input.dff(odIdx, ornIdx);    % response data block
    
    [expList, ia, ic] = unique(input.expID(odIdx), 'stable');   % find out how many trials each odor
    
    dffPoolMat  = reshape(dffPool,  [], length(expList));
    concPoolMat = reshape(concPool, [], length(expList)); 
    
    expIdx = find(~isnan(dffPoolMat(1, :)));
    
    concT123(:, 1:length(expIdx), i) = concPoolMat(:, expIdx);
    rspT123 (:, 1:length(expIdx), i) = dffPoolMat(:, expIdx);

end

%% prepare the initial guess
load(fullfile('results', 'initialParameters.mat'));
% 
for i = 1 : pairNum
    AVec0(i,1)  = aMatrix(sRow(i), sCol(i));
    x0Vec(i,1) = cMatrix(sRow(i), sCol(i));
end
% 
n = slop;
% nStd = .5;
% noiseStd = .4;

nStd = 1;
noiseStd = 1;

A = mean(AVec0);
AStd = std(AVec0);

% x0StdVec = .3;

x0StdVec = 1*ones(pairNum, 1);

%% fit
for iter = 1:5
    %first fit the global parameters
[A(iter+1),AStd(iter+1),n(iter+1),nStd(iter+1),noiseStd(iter+1),functionVal(iter+1),OutputInfo]...
    = FitMLE_global(rspT123, log10(concT123),...
    abs(normrnd(A(iter),AStd(iter)/2)),abs(normrnd(AStd(iter),AStd(iter)/2)),...
    x0Vec(:,iter),x0StdVec(:,iter),abs(normrnd(n(iter),nStd(iter))),abs(normrnd(nStd(iter),nStd(iter))),noiseStd(iter));
fprintf('Iter %d, A = %0.3f +- %0.3f, n = %0.3f +- %0.3f, function = %0.3f\n', iter, A(iter+1),...
    AStd(iter+1), n(iter+1),nStd(iter+1),functionVal(iter+1));
%then fit local parameters
[x0Vec(:,iter+1),x0StdVec(:,iter+1),funval(iter+1)] = ...
    FitMLE_local(rspT123,log10(concT123),A(iter+1),AStd(iter+1),x0Vec(:,iter),...
    x0StdVec(:,iter),n(iter+1),nStd(iter+1),noiseStd(iter+1));
fprintf('x0''s done, function = %0.3f\n',funval(iter+1))

end

%% save data and plot figures
% data is saved at ./results/MLEFit.mat.

