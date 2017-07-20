function [ dataMean, dataSEM, odorList, concList, infoORNList ] = GetAveNormData( )
%% Get the averaged matrix from the raw data

% read raw data
[ dataRaw, infoOdorRaw, ~, infoConcRaw, infoORNList] = ReadRawData();
% [~, dataRaw, infoOdorRaw, ~, infoConcRaw, infoORNList] = ReadRawData();

[odorList, ~] = unique(infoOdorRaw, 'stable');  %list of test odors
[concList, ~] = unique(infoConcRaw);            %list of test concentrations

% normalize data using the max actived neuron as the denormator.
dataNorm = dataRaw; %define the normlaized data
dataMean = zeros(length(odorList), length(infoORNList), length(concList));
dataSEM  = zeros(length(odorList), length(infoORNList), length(concList)); %standard error of the mean

for i = 1: length(odorList)
   od = odorList{i};
   indexTemp = strfind(infoOdorRaw, od);
   index = find(~cellfun('isempty', indexTemp));
   numTrial = length(index)/length(concList);
   if mod(numTrial, 1) ~= 0 % check if numTrial is an integer
       error(['missing data for odor ', od]);
   end
   
   for j = 1 : numTrial
        rowNum = index((j-1)*length(concList)+1 : j*length(concList));
        denom  = max(max(dataRaw(rowNum, :)));
        dataNorm(rowNum, :) = dataRaw(rowNum, :)./denom;
   end
   
   %average
   for k = 1 : length(concList)
       rowList = index(k : length(concList): (numTrial - 1)*length(concList) + k);
       dataMean(i, :, k) = mean(dataNorm(rowList, :));
       dataSEM(i, :, k)  = std(dataNorm(rowList, :))/sqrt(numTrial);
   end

end

end