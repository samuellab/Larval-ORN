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

%% go though the paris find all the non-active paris. 
fitMask = 99 * ones(length(input.Odor), length(input.ORN));

rMatSum = sum(rspTs, 3); 

fitMask(find(rMatSum == 0)) = 0; 
countPairZeros = length(find(rMatSum(:) == 0));
disp(['Count of non activated odor-ORN pairs:', num2str(countPairZeros), ...
    ', which is ',num2str(100*countPairZeros/numel(rMatSum)), ...
    '% of the total experiments.']);

%% pre-fit, look for pairs could be fitted well
% set up the fitting
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -11]; % setup the range and initial value of the variable
opts.Upper = [max(rspTs(:))*1.5 10 0];
opts.StartPoint = [4 4 -7];

[rowTotal, colTotal, ~] = size(rspTs);
cMatrix = NaN(rowTotal, colTotal);     % the thresholds for each odor_ORN pair, 'c' in the equation
aMatrix = NaN(rowTotal, colTotal);     % the saturated amplitude for each odor-ORN pair
r2Matrix= NaN(rowTotal, colTotal);     % the R-square value of the fitting

%%
disp('----------PRE-FIT----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(fitMask);
gfX = [];  gfY = [];    gfRC = [];  gfCoeff = [];	gfR2 = [];

for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft, opts);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    flag = coeff(1)>2 && coeff(3) < -5.8;
    if rSq > 0.9 && flag
       fitMask(ftRow(i), ftCol(i)) = 1;
       gfX = [gfX; log10(xx')];  gfY = [gfY; yy'];  gfRC = [gfRC; ftRow(i), ftCol(i)]; 
       gfCoeff = [gfCoeff; coeff];  gfR2 = [gfR2; rSq];
       
       fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            input.Odor{gfRC(end, 1)}, input.ORN{gfRC(end, 2)}, ...
            gfCoeff(end, 1), gfCoeff(end, 2), gfCoeff(end, 3),gfR2(end));
        
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);
        
        figure; plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{gfRC(end, 1)}, input.ORN{gfRC(end, 2)}]);
    end
end

disp(['Count of well-fitted saturated odor-ORN pairs:', num2str(length(gfR2))]);

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
plot(dataXEn(:), dataYEn(:), 'ok'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(xPlot, yPlot, 'r'); xlabel('Relative log10(c)'); ylabel('Norm.(\DeltaF/F)');
hold off;

for i = 1:length(ampVec)
    cMatrix(gfRC(i,1), gfRC(i,2)) = kdVec(i);
    aMatrix(gfRC(i,1), gfRC(i,2)) = ampVec(i);
    r2Matrix(gfRC(i,1), gfRC(i,2))= gfR2(i);
end

%%
ft2 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Bisquare';
opts2.Lower = [0 slop -11]; % setup the range and initial value of the variable
opts2.Upper = [max(rspTs(:))*1.5 slop 0];
opts2.StartPoint = [4 slop -7];


disp('----------FIT WITH KNOWN SLOP----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(fitMask == 99 );
gfX2 = [];  gfY2 = [];
gfRC2 = [];  gfCoeff2 = [];	gfR22 = [];
for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    flag = coeff(1)>2 && coeff(3) < -5;
    if rSq > 0.9 && flag
       fitMask(ftRow(i), ftCol(i)) = 2;
       gfX2 = [gfX2; log10(xx')];  gfY2 = [gfY2; yy'];  gfRC2 = [gfRC2; ftRow(i), ftCol(i)]; 
       gfCoeff2 = [gfCoeff2; coeff];  gfR22 = [gfR22; rSq];
       
       fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            input.Odor{gfRC2(end, 1)}, input.ORN{gfRC2(end, 2)}, ...
            gfCoeff2(end, 1), gfCoeff2(end, 2), gfCoeff2(end, 3),gfR22(end));
        
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);
        
        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{gfRC2(end, 1)}, input.ORN{gfRC2(end, 2)}]);
        
        cMatrix(ftRow(i), ftCol(i)) = coeff(3);
        aMatrix(ftRow(i), ftCol(i)) = coeff(1);
        r2Matrix(ftRow(i), ftCol(i))= rSq;
    end
end

disp(['Count of fitted odor-ORN pairs w/ predefined slop w/ EC50 smaller than -5:', ...
    num2str(length(gfR22))]);


%% take a look at the paris amplitude smaller than the setted threshold
[ftRow, ftCol] = find(fitMask);

disp('----------FIT WITH KNOWN SLOP,CHECK FAILED FITTINTS----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

for i = 1:length(ftRow)
    xx = squeeze(concTs(ftRow(i), ftCol(i), :));
    yy = squeeze(rspTs (ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    flag = coeff(1) <= 2 && coeff(1) > 0.4 && coeff(3) < -5;
    if rSq > 0.9 && flag

       fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', input.Odor{ftRow(i)}, ...
           input.ORN{ftCol(i)}, coeff(1), coeff(2), coeff(3),rSq);
        
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);
        
        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([input.Odor{ftRow(i)}, input.ORN{ftCol(i)}]);
        
    end
end






