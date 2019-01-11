% close all; clear; clc;
% warning('off','all');
% diary(fullfile('.', 'AnalysisResults', 'fitsingleORNdata_log.txt')); 
% diary on;

%% input file info
fileName = fullfile('data', 'Supplementary Table 1.csv');

% load Excel file
dataT = readtable(fileName);

[input.ORNs, ~, icORN] = unique(dataT.ORN, 'stable');     %ORN name

for i = 1:length(input.ORNs)
    list1 = find(icORN == i); 
    
    %odor name
    odorVec = unique(dataT.Odor(list1),'stable');       
    input.odors(i, :) = odorVec';   
    
    %concentration, df/f, exp_ID 
    for j = 1:length(input.odors(i,:))
        list2 = find(strcmp(input.odors{i, j}, dataT.Odor)) ;
        list12 = intersect(list1, list2);
        
        % concentration
        concVec = unique(dataT.Concentration(list12), 'stable'); input.concList{i, j} = concVec';
        
        % df/f
        input.dff{i, j} = transpose(reshape(dataT.DF_F(list12), length(input.concList{i, j}), []));
        
        % exp ID
        input.expID{i, j} = unique(dataT.Exp_ID(list12), 'stable'); 
    end
end



%% pool all the saturated data and format to 
 
concT = NaN(8, 12, 4); % concentration, odor-ORN pari, trial
rspT  = NaN(8, 12, 4);

for i = 1 : 2
    for j = 1:2
        
        idx = (i-1)*2 + j;
        
        cTemp = input.concList{i,j};
        rTemp = input.dff{i, j};
        
        trialNum = size(rTemp, 1);
        
        for k = 1:trialNum
            
            concT(:, k, idx) = cTemp';
            rspT(:, k, idx) = rTemp(k, :)';
        end
    end
end

%% prepare the initial guess
load(fullfile('AnalysisResults', 'fitSingleORNdataResults.mat'));

nTemp = zeros(4, 1);
for i = 1 : 2
    for j = 1 : 2
        idx = (i-1)*2 + j;
  
        
        AVec0(idx,1)  = mean(results.fitCoeffIdv{i, j}(:,1));
        x0Vec0(idx,1) = mean(results.fitCoeffIdv{i, j}(:,3));
        
        nTemp(idx) = results.fitCoeffEns{i, j}(2);
    end
end

AStdVec0 = ones(4, 1);
x0StdVec0 = 0.5 * ones(4, 1);

n0 = mean(nTemp);
nStd0 = 1;
noiseStd0 = 0.5;

%% fit

[AVec,AStdVec,x0Vec,x0StdVec,n,nStd,noiseStd, fval] = FitMLE_Ensemble(rspT, log10(concT), AVec0,AStdVec0,x0Vec0,x0StdVec0,n0,nStd0,noiseStd0)





%%
concb10 = log10(concT);
for ii = 1:size(rspT,3)
    conctrial = concb10(:,:,ii);
    cfit = linspace(min(conctrial(:)),max(conctrial(:)),20);
figure;
plot(squeeze(log10(concT(:,:,ii))),squeeze(rspT(:,:,ii)),'x-')
hold on;
plot(cfit,AVec(ii)./(1+exp(-n*(cfit-x0Vec(ii)))),'k','LineWidth',3)
plot(cfit,(AVec(ii)+AStdVec(ii))./(1+exp(-(n-nStd)*(cfit-(x0Vec(ii)-x0StdVec(ii))))),'k--')
plot(cfit,(AVec(ii)-AStdVec(ii))./(1+exp(-(n+nStd)*(cfit-(x0Vec(ii)+x0StdVec(ii))))),'k--')
end




