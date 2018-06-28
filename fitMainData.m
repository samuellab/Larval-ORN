fileName = fullfile('data', 'Supplementary Table 3.csv');

% load Excel file
dataT = readtable(fileName);

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');     %ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% 
