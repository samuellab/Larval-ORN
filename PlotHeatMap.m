clear; clc;

addpath(genpath(pwd));
%%
% diary(fullfile('.', 'AnalysisResults', 'PlotHeatMap.txt')); diary on;

% %% Read the .xlsx raw data file
% rawDataFile = fullfile('.', 'data', 'Supplementary Table 1.xlsx'); %define the .xlsx file
% sheet = 1;
% xlRange = 'A1:X611';
% [~,~,raw] = xlsread(rawDataFile, sheet, xlRange);
% 
% %define the experiment information list
% infoOdorRaw  = raw(2:end, 1);
% infoExpIDRaw = raw(2:end, 2);
% infoConcRaw  = cell2mat(raw(2:end, 3));
% 
% infoORNListRaw  = raw(1, 4:end);                   %list of idendtified neurons
% [odorList, ~] = unique(infoOdorRaw, 'stable');  %list of test odors
% [concList, ~] = unique(infoConcRaw);            %list of test concentrations
% 
% % extract the data matrix
% dataRaw = cell2mat(raw(2:end, 4:end));
% 
% % find out the unidentified neurons (have no responses)
% index0 = sum(dataRaw, 1)==0;
% 
% % get rid of unidentified neurons in the data and ORN list
% index1 = sum(dataRaw, 1)~=0;
% infoORNList = infoORNListRaw(index1);
% data = dataRaw(:, index1);
% 
% disp('----------Information About the Raw Dose Response Data----------');
% fprintf('%30s %.0f \n', 'Number of odor:', length(odorList));
% fprintf('%30s %.0f \n', 'Concentrations each odor:', length(concList));
% fprintf('%30s %.0f \n', 'Number of larvae imaged:', (length(raw)-1)/length(concList));
% fprintf('%30s %.1f \n', 'Average trials each condition:', length(raw)/(length(odorList)*length(concList)));
% fprintf('%30s %.0f \n', 'Number of responded ONRs:', length(infoORNList));
% fprintf('%30s %5s %5s %5s \n', 'ONRs no responses:',  infoORNListRaw{index0});
% 
% %% average data cross trials
% dataRawAve = zeros(length(odorList), length(infoORNList), length(concList));
% dataRawSEM = zeros(length(odorList), length(infoORNList), length(concList));
% 
% for i = 1: length(odorList)
%    od = odorList{i};
%    indexTemp = strfind(infoOdorRaw, od);
%    index = find(~cellfun('isempty', indexTemp));
%    numTrial = length(index)/length(concList);
%    if mod(numTrial, 1) ~= 0 % check if numTrial is an integer
%        error(['missing data for odor ', od]);
%    end
% 
%    for j = 1 : length(concList)
%        rows =  index(j: length(concList): (numTrial-1)*length(concList)+j);
%        for k = 1:length(infoORNList)
%            dataRawAve(i, k, j) = mean(data(rows, k));
%            dataRawSEM(i, k, j) = std(data(rows, k))/sqrt(numTrial);
%        end
%    end
% end
% 
% % save data to file
% save(fullfile('.', 'data', 'AveRawDataMatrix.mat'), 'dataRawAve', ...
%     'dataRawSEM', 'odorList', 'concList', 'infoORNList');

%% Visualize the data, plot the heatmap
% in the paper, we decide to use 18 odors instead of 19. So we removed one.
load(fullfile('.', 'data', 'AveRawDataMatrix2ndRound.mat'));

odorListV = odorList;
dataRawAveV = dataMean; 

% odorListV(5) = [];  dataRawAveV(5, :, :)= [];

% Simulated annealing algorithm is used to search a combination of odor and
% ORN order to make the matrix looks 'organized'.
reproduceFig2 = 0;
if reproduceFig2
    % Because SA is a stochastic algorithm, the result may be different for each 
    % run. To reproduce the figure in the paper, use the following orders: 
%     odorOrder = [17 12 15 2 10 3 4 16 9 1 18 8 6 5 7 13 14 11]; 
%     ORNOrder  = [16 17 5 2 14 12 11 1 4 7 10 13 15 9 18 8 6 3];
    odorOrder = [19,33,12,32,27,15,14,7,8,9,6,22,24,31,30,29,1,5,10,17,28,3,4,23,20,26,13,25,16,2,11,21,34,18]; 
    ORNOrder  = [19,21,3,6,8,14,10,16,9,7,11,4,1,13,12,2,20,5,15,17,18];
else
    % Use simulated annealing algorithm to search a order of row and col
    [odorOrder, ORNOrder] = GetBestOrder( dataRawAveV );
end

% organize the matrix using the odor and ORN order
mNewOrderTemp = dataRawAveV;
for i = 1:length(odorOrder)
    mNewOrderTemp(i,:,:) = dataRawAveV(odorOrder(i),:,:);
end
mNewOrder = mNewOrderTemp;
for i = 1:length(ORNOrder)
    mNewOrder(:,i,:) = mNewOrderTemp(:,ORNOrder(i), :);
end

%% draw the heatmap
for i = 1 : length(concList)
    data2D = mNewOrder(:,:,i);
    figure;
    set(gcf, 'Position', [100 250 560 700] );
    a = axes;
    set(a, 'CLim', [min(min(min(mNewOrder))) 1]);
    imagesc(data2D); 
    set(gca,'XTick',1:length(ORNOrder));
    set(gca,'XTickLabel',ORNList(ORNOrder));
    set(gca,'xaxisLocation','top');
    set(gca,'YTick',1:length(odorOrder));
    set(gca,'YTickLabel',odorListV(odorOrder));
    ax = gca; ax.XTickLabelRotation = 45;
    title(num2str(concList(i)));
    colormap(jet);
    caxis([min(min(min(mNewOrder))) max(max(max(mNewOrder)))]);
end

% diary off;