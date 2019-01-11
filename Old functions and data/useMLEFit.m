fileName = fullfile('data', 'Supplementary Table 3.csv');

dataT = readtable(fileName);    % load Excel file

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');	%ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% mark the type of each curve
maskPlot = SetCurveType(34, 21); % manually define the type of the dose-response curves

%% pool all the saturated data and format to 
[sRow, sCol] = find(maskPlot == 1);

pairNum = length(sRow);

concT1 = NaN(5, 10, pairNum);
rspT1  = NaN(5, 10, pairNum);

for i = 1 : pairNum
    odIdx  = find(input.indOdor == sRow(i)); % row index of each odor
    ornIdx = sCol(i);
    
    concPool = input.concList(odIdx);       % concentration list
    dffPool  = input.dff(odIdx, ornIdx);    % response data block
    
    [expList, ia, ic] = unique(input.expID(odIdx), 'stable');   % find out how many trials each odor
    
    dffPoolMat  = reshape(dffPool,  [], length(expList));
    concPoolMat = reshape(concPool, [], length(expList)); 
    
    expIdx = find(~isnan(dffPoolMat(1, :)));
    
    concT1(:, 1:length(expIdx), i) = concPoolMat(:, expIdx);
    rspT1 (:, 1:length(expIdx), i) = dffPoolMat(:, expIdx);

end

%% prepare the initial guess
load('fitResults.mat');

for i = 1 : pairNum
    AVec0(i,1)  = aMatrix(sRow(i), sCol(i));
    x0Vec0(i,1) = cMatrix(sRow(i), sCol(i));
end

AStdVec0 = ones(pairNum, 1);
x0StdVec0 = 0.5 * ones(pairNum, 1);

% n0 = slop;
n0 = 4.8;
nStd0 = 0.5;
noiseStd0 = 1.8;

%% fit

[AVec,AStdVec,x0Vec,x0StdVec,n,nStd,noiseStd] = FitMLE_Ensemble(rspT1, concT1, AVec0,AStdVec0,x0Vec0,x0StdVec0,n0,nStd0,noiseStd0)















