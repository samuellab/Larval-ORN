addpath(fullfile('.', 'tools'));

%%
fileName = fullfile('data', 'Supplementary Table 3.csv');

dataT = readtable(fileName);    % load Excel file

[input.Odor, ~, input.indOdor] = unique(dataT.Odor, 'stable');	%ORN name
input.ORN = dataT.Properties.VariableNames(4:end);
input.concList = dataT.Concentration;
input.expID = dataT.Exp_ID;
input.dff = table2array(dataT(:, 4:end));

%% average data cross trials for each odor-ORN pair
concTs = zeros(length(input.Odor), length(input.ORN), 5);
rspTs  = concTs;

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
            concTs(i, j, :) = concBlock(:, 1);
            rspTs(i, j, :) = mean(dffBlock, 2);
        else
            error('Trials do not share the same set of odor concentration.');
        end
    end
end

%% define the curve type
[rowTotal, colTotal, ~] = size(concTs);
maskPlot = SetCurveType(rowTotal, colTotal); % manually define the type of the dose-response curves

%% pre-fit, see which pair could be fitted well
% set up the fitting
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -11]; % setup the range and initial value of the variable
opts.Upper = [max(rspTs(:))*1.5 10 0];
opts.StartPoint = [4 4 -7];

cMatrix = NaN(rowTotal, colTotal);     % the thresholds for each odor_ORN pair, 'c' in the equation
aMatrix = NaN(rowTotal, colTotal);     % the saturated amplitude for each odor-ORN pair
r2Matrix= NaN(rowTotal, colTotal);     % the R-square value of the fitting

%% fit each saturate pairs
disp('----------PRE-FIT----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

plotFlag = 0;

[ftRow, ftCol] = find(maskPlot == 1);
gfX = [];  gfY = [];  gfRC = [];  gfCoeff = [];	gfR2 = [];
for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft, opts);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);

    gfX = [gfX; log10(xx')];  gfY = [gfY; yy'];  gfRC = [gfRC; ftRow(i), ftCol(i)]; 
    gfCoeff = [gfCoeff; coeff];  gfR2 = [gfR2; rSq];

    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', input.Odor{gfRC(end, 1)}, ...
        input.ORN{gfRC(end, 2)},coeff(1), coeff(2), coeff(3), rSq);

    if plotFlag
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{gfRC(end, 1)}, input.ORN{gfRC(end, 2)}]);
    end
end

%% ensemble fit of all these good fitted curves
slop0 = median(gfCoeff(:, 2)); ampVec0 = gfCoeff(:, 1); kdVec0 = gfCoeff(:, 3);

disp('----------FIT CURVE ENSEMBLE:----------');
fprintf('%-5s\t%-5s\t\n', 'Slop', 'R^2');
[slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(gfX, gfY, hillEq, slop0, ampVec0, kdVec0);

fprintf('%.2f\t%.2f\t\n', slop, rSquare);

% plot the data 
dataXEn = gfX -  repmat(kdVec,  1, length(gfX(1,:)));
dataYEn = gfY ./ repmat(ampVec, 1, length(gfY(1,:)));

% save into results
results.fitCoeffFMS = [ampVec, repmat(slop, length(ampVec), 1), kdVec];

figure; 
plot(10.^(dataXEn'), dataYEn', 'o'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(10.^xPlot, yPlot, 'r'); xlabel('Relative Concentration'); ylabel('Norm.(\DeltaF/F)');
set(gca, 'XScale', 'log')
hold off;
for i = 1:length(ampVec)
    cMatrix(gfRC(i,1), gfRC(i,2)) = kdVec(i);
    aMatrix(gfRC(i,1), gfRC(i,2)) = ampVec(i);
    r2Matrix(gfRC(i,1), gfRC(i,2))= gfR2(i);
end

%% fit with fixed slop
ft2 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Bisquare';
opts2.Lower = [0 slop -11]; % setup the range and initial value of the variable
opts2.Upper = [max(rspTs(:))*1.5 slop 0];
opts2.StartPoint = [4 slop -7];

disp('----------FIT WITH KNOWN SLOP----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(maskPlot == 2);
% gfX2 = [];  gfY2 = [];
% gfRC2 = [];  gfCoeff2 = [];	gfR22 = [];

for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    aMatrix(ftRow(i), ftCol(i)) = coeff(1);
    cMatrix(ftRow(i), ftCol(i)) = coeff(3);
    r2Matrix(ftRow(i), ftCol(i))= rSq;
       
    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', input.Odor{ftRow(i)}, ...
        input.ORN{ftCol(i)}, coeff(1), coeff(2), coeff(3), rSq);
        
	if plotFlag == 1
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{ftRow(i)}, input.ORN{ftCol(i)}]);
    end
end

% %% show the distribution of the amplitudes of the type 3 curves.
% [ftRow, ftCol] = find(maskPlot == 3);
% ymax3 = [];
% for i = 1:length(ftRow)
%     yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
%     
%     ymax3 = [ymax3; max(yy(:))];
% end


%% Estimate the maximum values for each Odor

% find out all the non-nan values in the aMatrix and plot the histogram
% ampKnown=aMatrix(~isnan(aMatrix));
% histogram(ampKnown);

% meanA = nanmean(aMatrix(:));

maxAofOdor = zeros(rowNum, 1); % define parameter to store maximum amplitude 
lowerBound = 1; % set a lower bound of the maximum amplitude

disp('----------FIND AMPLITUDE OF EACH ODOR FOR TYPE 3 CURVE FITTING----------');
fprintf('%30s\t%-5s\t\n', 'Odor Name', 'Amplitude');

maxProj = max(rspTs, [], 3); %max of the data along the concentration
 
for i = 1:rowNum
    index = ~isnan(aMatrix(i, :)); %find the nun-zero elements of each row (odor), from the saturated dataset
    maxSeqTemp = aMatrix(i, index);
    maxSeq = maxSeqTemp(maxSeqTemp>lowerBound); %select the max values higher than the lower bound
        
    maxEst = mean(maxSeq); % use the mean as the maximum
%     if ~isnan(maxEst) && maxEst>lowerBound
        maxAofOdor(i) = maxEst;
%     else %if there is no qualified elements, find the maximum in the raw data
%         indexBackup = find(maxProj(i, :)>lowerBound);
%         seqBackup = maxProj(i, indexBackup);
%         maxEstBackup = mean(seqBackup);
%         if isnan(maxEstBackup)
%             maxAofOdor(i) = max(maxProj(i, :));
%         else
%             maxAofOdor(i) = maxEstBackup;
%         end
%     end
    
	fprintf('%30s\t%.2f\t\n', input.Odor{i}, maxAofOdor(i));
end

%% plot selected odor's dose-response map
selOdorIdx = 1;
figure; 
for i = 1:colNum
    xx = squeeze(concTs(selOdorIdx, i, :));
    yy = squeeze(rspTs(selOdorIdx, i, :));
    subplot(3, 7, i);
    plot(log10(xx), yy, 'o-'); title(input.ORN{i});
end


%% fit with fixed slop and amplitude
ft3 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Robust = 'Bisquare';
% opts3.Lower = [meanA slop -11]; % setup the range and initial value of the variable
% opts3.Upper = [meanA slop 0];
% opts3.StartPoint = [meanA slop -7];

disp('----------FIT WITH KNOWN SLOP----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(maskPlot == 3);

plotFlag = 0;
for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));

    if max(yy) > maxAofOdor(ftRow(i))
        [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);
        
    else
        % setup the range and initial value of the variable
        opts3.Lower = [maxAofOdor(ftRow(i)) slop -11]; 
        opts3.Upper = [maxAofOdor(ftRow(i)) slop 0];
        opts3.StartPoint = [maxAofOdor(ftRow(i)) slop -7];
        
        [fitresult, gof] = fit(log10(xx), yy, ft3, opts3);   %fit
    end
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    aMatrix(ftRow(i), ftCol(i)) = coeff(1);
    cMatrix(ftRow(i), ftCol(i)) = coeff(3);
    r2Matrix(ftRow(i), ftCol(i))= rSq;
       
    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', input.Odor{ftRow(i)}, ...
        input.ORN{ftCol(i)}, coeff(1), coeff(2), coeff(3), rSq);
        
	if plotFlag == 1
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{ftRow(i)}, input.ORN{ftCol(i)}]);
    end
end


% %%
% 
% % compare by row,
% varRow = [];
% for i = 1: rowNum
%     rowData = aMatrix(i, :);
%     rowData = rowData(~isnan(rowData));
%     
%     rowData = rowData - mean(rowData);
%     varRow = [varRow, rowData];
% end
% std(varRow)
% 
% 
% varCol = [];
% for j = 1: colNum
%     colData = aMatrix(:, j);
%     colData = colData(~isnan(colData));
%     
%     colData = colData - mean(colData);
%     varCol = [varCol; colData];
% end
% std(varCol)
% 

%% show the failed fittings
gofThld = 0.5;

disp('----------FAILED FITTINGS----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

for i = 1: rowTotal
    for j = 1: colTotal
        if r2Matrix(i, j) < gofThld
            fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', input.Odor{i}, ...
                input.ORN{j}, aMatrix(i,j), slop, cMatrix(i,j), r2Matrix(i,j));
            
            if plotFlag == 1
                xx = squeeze(concTs(i, j, :));
                yy = squeeze(rspTs(i, j, :));
    
                xP = linspace(min(log10(xx)), max(log10(xx)), 50);
                yP = hillEq( aMatrix(i,j), slop, cMatrix(i,j), xP);

                figure; 
                plot(log10(xx), yy, 'ok'); hold on;
                plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
                title([input.Odor{i}, input.ORN{j}]);
            end
        end
    end
end

%% show the EC50 matrix and amlitude matrix

%replace NaN in the EC50 matrix with 0
cMatrix(isnan(cMatrix)) = 0;

% apply the sequence of ORN and odor to order elements of the matrix 
% odorOrder = [17 12 15 2 10 3 4 16 9 1 18 8 6 5 7 13 14 11]; % the order is consistant to figure 2
% ORNOrder = [16 17 5 2 14 12 11 1 4 7 10 13 15 9 18 8 6 3];

% odorOrder = [19,33,12,32,27,15,14,7,8,9,6,22,24,31,30,29,1,5,10,17,28,3,4,23,20,26,13,25,16,2,11,21,34,18]; 
% ORNOrder  = [19,21,3,6,8,14,10,16,9,7,11,4,1,13,12,2,20,5,15,17,18];

odorOrder = [18,34,21,2,11,16,25,13,26,17,20,28,23,4,3,10,5,1,29,27,24,30,31,22,6,9,8,7,14,15,32,12,19,33]; 
ORNOrder  = [18,17,15,5,20,2,1,12,13,4,11,7,16,9,10,14,8,6,3,21,19];

[rowNum, colNum] = size(cMatrix);
newMStep1 = -cMatrix;
for i = 1:rowNum
    newMStep1(i,:) = -cMatrix(odorOrder(i),:);
end
newMStep2 = newMStep1;
for i = 1:colNum
    newMStep2(:,i) = newMStep1(:,ORNOrder(i));
end

% draw the heat map of the EC50 matrix
ec50Map = newMStep2;
figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 610, 420]);
imagesc(ec50Map); 
set(gcf, 'Position', [100 250 560 700])
set(gca, 'CLim', [0 max(ec50Map(:))]);
set(gca,'XTick',1:colNum);
set(gca,'XTickLabel',input.ORN(ORNOrder));
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:rowNum);
set(gca,'YTickLabel',input.Odor(odorOrder));
set(gca, 'XTickLabelRotation', 45);
cmp = colormap(jet); cmp(1,:) = [0 0 0];
colormap(cmp); c = colorbar; 
c.TickLabels{1} = 'NaN'; c.Label.String = '-log10(EC50)';
title('EC50'); 


% show the amplitude matrix
aMatShow = aMatrix;    aMatShow(isnan(aMatShow)) = 0;

newAStep1 = aMatShow;
for i = 1:rowNum
    newAStep1(i,:) = aMatShow(odorOrder(i),:);
end
newAStep2 = newAStep1;
for i = 1:colNum
    newAStep2(:,i) = newAStep1(:,ORNOrder(i));
end

% draw the heat map of the EC50 matrix
aMap = newAStep2;
figure;  pos = get(gcf, 'pos'); set(gcf, 'pos', [pos(1), pos(2), 610, 420]);
imagesc(aMap); 
set(gcf, 'Position', [100 250 560 700])
set(gca, 'CLim', [0 max(aMap(:))]);
set(gca,'XTick',1:colNum);
set(gca,'XTickLabel',input.ORN(ORNOrder));
set(gca,'xaxisLocation','top');
set(gca,'YTick',1:rowNum);
set(gca,'YTickLabel',input.Odor(odorOrder));
set(gca, 'XTickLabelRotation', 45);
cmp = colormap(jet); cmp(1,:) = [0 0 0];
colormap(cmp); c = colorbar; 
c.TickLabels{1} = 'NaN'; c.Label.String = '-log10(EC50)';
title('Amp'); 


%%
% save('fitResults.mat', 'aMatrix', 'cMatrix', 'r2Matrix', 'slop');





