function [ dataMean, dataSEM, odorList, concList, infoORNList ] = GetAveData( )
%% Get the averaged matrix from the raw data

% read raw data
[ dataRaw, infoOdorRaw, ~, infoConcRaw, infoORNList] = ReadRawData();
% [~, dataRaw, infoOdorRaw, ~, infoConcRaw, infoORNList] = ReadRawData();

[odorList, ~] = unique(infoOdorRaw, 'stable');  %list of test odors
[concList_full, ~] = unique(infoConcRaw);            %list of test concentrations

%
if length(concList_full) > 5
    concList = concList_full(end-4 : end);
else
    concList = concList_full;
end

%% mean and SEM of the data
dataMean = zeros(length(odorList), length(infoORNList), length(concList));
dataSEM  = zeros(length(odorList), length(infoORNList), length(concList)); %standard error of the mean

for i = 1: length(odorList) %go through odors
   od = odorList{i};
   indexTemp = strfind(infoOdorRaw, od);
   index_odor = find(~cellfun('isempty', indexTemp));
   
   % go through concentrations
   for j = 1 : length(concList)
       index_con = find(infoConcRaw(index_odor) == concList(j));
       rowList = index_odor(index_con);
       dataBlock = dataRaw(rowList, :);
       for k = 1:length(infoORNList) %go though ORNs
           dataVec = dataBlock(:, k);
           index_nnan = find(~isnan(dataVec));
           dataVecPure = dataVec(index_nnan);
           dataMean(i, k, j) = mean(dataVecPure);
           dataSEM(i, k, j)  = std(dataVecPure)/sqrt(length(dataVecPure));
       end
   end
end

% %% visualize the data, roughly
% for i = 1:length(concList)
%     dataMat = dataMean(:, :, i);
%     figure; imagesc(dataMat); title(['mean ', num2str(concList(i))]);
%     
%     semMat = dataSEM(:,:,i);
%     figure; imagesc(semMat); title(['SEM ', num2str(concList(i))]);
% end

end