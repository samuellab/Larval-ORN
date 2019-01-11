%% read raw data
% [rawT, dataRaw, infoOdorRaw, infoExpIDRaw, infoConcRaw, infoORNList] = ReadRawData();
load(fullfile('data', 'doseResponseData.mat'));

% dataRawTemp=dataRaw(1:530, :);

% [odorList, ~] = unique(infoOdorRaw, 'stable');  %list of test odors
% [concList, ~] = unique(infoConcRaw);            %list of test concentrations
concList = concHm;
dataRaw = dffHm;

% dataNorm = dataRaw; %define the normlaized data
% for i = 1: length(odorList)
%    od = odorList{i};
%    indexTemp = strfind(infoOdorRaw, od);
%    index = find(~cellfun('isempty', indexTemp));
%    numTrial = length(index)/length(concList);
%    if mod(numTrial, 1) ~= 0 % check if numTrial is an integer
%        error(['missing data for odor ', od]);
%    end
% 
%    for j = 1 : numTrial
%         rowNum = index((j-1)*length(concList)+1 : j*length(concList));
%         denom  = max(max(dataRaw(rowNum, :)));
%         dataNorm(rowNum, :) = dataRaw(rowNum, :)./denom;
%    end
% end

classTemp = rawT{:, 1};

%% 
numTop = 3;
dataSelect = [];
class= cell(length(dataRaw)*numTop/5, 1);
for i = 1:length(dataRaw)/5
    indexL = (i-1)*5 + 5-numTop+1;
    indexR = (i-1)*5 + 5;
    dataSelect = [dataSelect; dataRaw(indexL:indexR, :)];
    class((i-1)*numTop+1: (i-1)*numTop+numTop) = classTemp(indexL:indexR);
end

labels =  unique(infoOdorRaw);
ORNs = dataSelect;

%%
totalRun = 1;
confMatAll = zeros(length(odorList), length(odorList), totalRun);
ORNdata = ORNs;
for i=1:totalRun
    i
    cv = cvpartition(class,'holdout',0.2);

    ORNdata_train = ORNdata(training(cv),:);
    class_train = class(training(cv));

    ORNdata_test = ORNdata(test(cv),:);
    class_test = class(test(cv));

    %% train linear discriminant analysis classifier
    mdl = ClassificationDiscriminant.fit(ORNdata_train, class_train);

    %% predict classes using the LDA model
    predicted_classes = predict(mdl,ORNdata_test);

    %% compute and visualize confusion matrix
    conf_mat = confusionmat(class_test,predicted_classes);
    confMatAll(:,:,i) = conf_mat;
end

confMatAve = sum(confMatAll, 3);
accr = trace(confMatAve)/sum(sum(confMatAve));
confMatAveNor = confMatAve;

for i =1:length(odorList)
    confMatAveNor(i,:)= confMatAve(i,:)./sum(confMatAve(i,:));
end

figure;
heatmap(confMatAveNor, labels, labels, 0,   'TickAngle', 45, 'Colormap', 'red',...
             'ShowAllTicks',1,  'Colorbar', true);
title(['Top ', num2str(numTop), ' Concentration, A=', num2str(accr)]);
disp(['Accuracy: ', num2str(accr)]);

%%
deleteTimes = 17;
repTimes = 20;

accrMat = zeros(deleteTimes, repTimes);
for jj = 1:deleteTimes
    jj
    for kk = 1:repTimes
        colNum = length(odorList) - jj;
        rIndex = randperm(length(infoORNList), colNum);
        ORNdata = ORNs(:, rIndex);

        %% partition data into training and testing sets
        %c = cvpartition(group,'HoldOut',p) randomly partitions observations into a
        %training set and a test set with stratification, using the class information
        %in group; that is, both training and test sets have roughly the same class 
        %proportions as in group.
        totalRun = 2;
        confMatAll = zeros(length(odorList), length(odorList), totalRun);
        for i=1:totalRun
             cv = cvpartition(class,'holdout',0.2);

            ORNdata_train = ORNdata(training(cv),:);
            class_train = class(training(cv));

            ORNdata_test = ORNdata(test(cv),:);
            class_test = class(test(cv));

            %% train linear discriminant analysis classifier
            mdl = ClassificationDiscriminant.fit(ORNdata_train, class_train);

            %% predict classes using the LDA model
            predicted_classes = predict(mdl,ORNdata_test);

            %% compute and visualize confusion matrix
            conf_mat = confusionmat(class_test,predicted_classes);
            confMatAll(:,:,i) = conf_mat;
        end
        confMatAve = sum(confMatAll, 3);
        accr = trace(confMatAve)/sum(sum(confMatAve));
        accrMat(jj, kk) = accr;
        
%         figure;
%         heatmap(confMatAve, labels, labels, 1, 'TickAngle', 45, 'Colormap', 'red',...
%              'ShowAllTicks',1,'UseLogColorMap', true, 'Colorbar', true);
%         title(['Top ', num2str(numTop), ' Concentration, A=', num2str(accr)]);
%         disp(['Accuracy: ', num2str(accr)]);

        
    end
end
aveAccrMat = mean(accrMat, 2);
figure;
% plot(1:5,accrList(end:-1:1), 'o'); xlabel('Concentration Level'); ylabel('accuracy');
for i =1:repTimes
    plot(1:deleteTimes, accrMat(:, i), 'o'); hold on;
end
plot(1:deleteTimes, aveAccrMat, '*k'); hold on;
hold off; xlabel('# of Removed ORNs'); ylabel('Accuracy');

