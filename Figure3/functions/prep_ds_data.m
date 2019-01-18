function [t123idx,sRow,sCol,concT123,rspT123] = prep_ds_data()
%% Load the raw dose response data and initial guesses from least squares fitting
addpath(genpath(pwd));

fileName =  which('Data S1.csv');

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

end