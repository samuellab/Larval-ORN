function [ dataRaw, infoOdorRaw, infoExpIDRaw, infoConcRaw, infoORNList] = ReadRawData()
%% Read the .xlsx data file in the subfolder of 'data',
% no input varialbes,
% outputs are the raw data matrix and information of the experiments

filename = fullfile('..', 'data_2nd_round', 'Supplementary Table 1_plusNewData.xlsx'); %define the .xlsx file
sheet = 1;
xlRange = 'A1:X1161';
[~,~,raw] = xlsread(filename,sheet,xlRange);

%define the experiment information list
infoOdorRaw  = raw(2:end, 1);
infoExpIDRaw = raw(2:end, 2);
infoConcRaw  = cell2mat(raw(2:end, 3));
infoORNList  = raw(1, 4:end);

%%
%extract the data matrix
A = importdata(filename);
B=A.data.Sheet1;
dataRaw = B (:, 3:end);

%% clean the data
% % find out the unidentified (no activity) neurons
% index0 = sum(dataRaw, 1)==0;
% disp('List of unidentified neurons:');
% disp(infoORNList(index0));

% remove the unidentified neurons in the data and ORN list, and rewrite the varaibles
% index1 = sum(dataRaw, 1)~=0;
% infoORNList = infoORNList(index1);
% dataRaw = dataRaw(:, index1);

% %% organize the informaiton and data into a table
% tColHeader = [{'Odor', 'Exp_ID', 'Concentration'}, infoORNList];
% tData = cell(length(dataRaw), length(tColHeader));
% tData(:, 1)     = infoOdorRaw;
% tData(:, 2)     = infoExpIDRaw;
% tData(:, 3)     = num2cell(infoConcRaw);
% tData(:, 4:end) = num2cell(dataRaw);
% 
% rawT = cell2table( tData , 'VariableNames',tColHeader); %define the table
% 
% %%
% infoORNList{1}  = 'Or33b-47a';
% infoORNList{18} = 'Or94a-94b';

end