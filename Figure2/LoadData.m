%% read csv file
fileName = fullfile('data', 'Data S1.csv');

rawDataTable = readtable(fileName);    % load Excel file

[input.Odor, ~, input.indOdor] = unique(rawDataTable.Odor, 'stable');	%ORN name
input.ORN = rawDataTable.Properties.VariableNames(4:end);
input.concList = rawDataTable.Concentration;
input.expID = rawDataTable.Exp_ID;
input.dff = table2array(rawDataTable(:, 4:end));

%% average data cross trials for each odor-ORN pair
conc = zeros(length(input.Odor), length(input.ORN), 5); dff = conc; dffSEM = conc;

odorList = input.Odor;  ORNList = input.ORN;

% for a better display in Matlab plot
ORNList = cellfun(@(x) strrep(x, '_', '/'), ORNList, 'UniformOutput', false);

for i = 1 : length(input.Odor)
    odI = find(input.indOdor == i); % row index of each odor
    
    concPool = input.concList(odI); % concentration list
    dffPool = input.dff(odI, :);    % response data block
    
    [expList, ia, ic] = unique(input.expID(odI), 'stable');   % find out how many trials each odor
    
    dffPoolTensor = reshape(dffPool', length(input.ORN), [], length(expList));
    concPoolMat = reshape(concPool, [], length(expList)); 
    
    for j = 1 : length(input.ORN)
        dffMat = squeeze(dffPoolTensor(j, :, :));
        
        % each colume is a trial, find trilas are not NaN
        colIdx = find(~isnan(dffMat(1, :)));
        
        % check if the concentration list is the same for these trials
        concBlock =concPoolMat(:, colIdx);
        dffBlock = dffMat(:, colIdx);
        if isequal(concBlock(:, 1:end-1), concBlock(:, 2:end))
            conc(i, j, :) = concBlock(:, 1);
            dff(i, j, :) = mean(dffBlock, 2);
            dffSEM(i, j, :) = std(dffBlock, 0, 2)/sqrt(size(dffBlock, 2));
        else
            error('Trials do not share the same set of odor concentration.');
        end
    end
end

%% for ultral-sensitive pairs, fill the missing data with highest concentration data
odorCount = size(conc, 1);      ORNCount  = size(conc, 2);

mark = sum(conc - repmat(conc(1,1,:), [odorCount, ORNCount, 1]), 3);

[rowFill, colFill] = find(mark);

r = randi(odorCount);	c = randi(ORNCount);    concTarget = conc(r, c, :);
while mark(r, c) ~= 0
    r = randi(odorCount); c = randi(ORNCount);
    concTarget = conc(r, c, :);
end

concHm = squeeze(concTarget); % 'Hm' stands for Heat map
dffHm = dff;  dffSEMHm = dffSEM;    

for i = 1 : length(rowFill)
    concCrt = conc(rowFill(i), colFill(i), :);
    dffCrt  = dff(rowFill(i), colFill(i), :);
    dffSEMCrt = dffSEM(rowFill(i), colFill(i), :);
    
    [C, ia, ib] = intersect(concCrt, concTarget);

    dffTgt = repmat(dffCrt(end), size(dffCrt));  dffTgt(ib) = dffCrt(ia);
    dffSEMTgt = repmat(dffSEMCrt(end), size(dffCrt)); dffSEMTgt(ib) = dffSEMCrt(ia);
    
    dffHm(rowFill(i), colFill(i), :) = dffTgt;
    dffSEMHm(rowFill(i), colFill(i), :) = dffSEMTgt;
end

%% correct typos on odor name
odorList{find(strcmp(odorList,'4-methylcyclohexane' ))} = '4-methylcyclohexanol ';
odorList{find(strcmp(odorList,'4-pheny-2-butanol' ))} = '4-phenyl-2-butanol';

%% save data
save(fullfile('results', 'doseResponseData.mat'), 'rawDataTable', ...
	'odorList', 'ORNList', 'conc', 'dff', 'dffSEM', 'concHm', 'dffHm', 'dffSEMHm');