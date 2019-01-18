addpath(genpath(pwd));

% load data
load(fullfile('data', 'doseResponseData.mat'));

% define the curve type
[rowTotal, colTotal, ~] = size(conc);
maskPlot = SetCurveType(rowTotal, colTotal); % manually define the type of the dose-response curves

% set up the fitting
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));

ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [0 0 -11]; % setup the range and initial value of the variable
opts.Upper = [max(dff(:))*1.5 10 0];
opts.StartPoint = [4 4 -7];

cMatrix = NaN(rowTotal, colTotal);     % the thresholds for each odor_ORN pair, 'c' in the equation
aMatrix = NaN(rowTotal, colTotal);     % the saturated amplitude for each odor-ORN pair
r2Matrix= NaN(rowTotal, colTotal);     % the R-square value of the fitting
stdMatrix=NaN(rowTotal, colTotal);

%% fit each saturated curve
disp('----------PRE-FIT----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

plotFlag = 0; % flg to plot each individual curve

[ftRow, ftCol] = find(maskPlot == 1); % maskPlot == 1 means, the curve is saturated
gfX = [];  gfY = [];  gfRC = [];  gfCoeff = [];	gfR2 = [];
for i = 1:length(ftRow)
    xx = squeeze(conc(ftRow(i), ftCol(i), :));
    yy = squeeze(dff(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft, opts);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);

    gfX = [gfX; log10(xx')];  gfY = [gfY; yy'];  gfRC = [gfRC; ftRow(i), ftCol(i)]; 
    gfCoeff = [gfCoeff; coeff];  gfR2 = [gfR2; rSq];

    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', odorList{gfRC(end, 1)}, ...
        ORNList{gfRC(end, 2)},coeff(1), coeff(2), coeff(3), rSq);

    if plotFlag
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([odorList{gfRC(end, 1)}, ORNList{gfRC(end, 2)}]);
    end
end

%% fit the saturated curves together. Constrain they share the same optimized slop. Plot Figure 3A and Figure S4A
slop0 = median(gfCoeff(:, 2)); ampVec0 = gfCoeff(:, 1); kdVec0 = gfCoeff(:, 3);

disp('----------SEARCH BEST SLOP FOR ALL SATURATED CURVES:----------');
fprintf('%-5s\t%-5s\t\n', 'Slop', 'R^2');
[slop, ampVec, kdVec, rSquare] = EnsembleMiniSearch(gfX, gfY, hillEq, slop0, ampVec0, kdVec0);

fprintf('%.2f\t%.2f\t\n', slop, rSquare);

% plot the data 
dataXEn = gfX -  repmat(kdVec,  1, length(gfX(1,:)));
dataYEn = gfY ./ repmat(ampVec, 1, length(gfY(1,:)));

% save into results
results.fitCoeffFMS = [ampVec, repmat(slop, length(ampVec), 1), kdVec];

% plot Figure 3A.
figure; 
plot(10.^(dataXEn'), dataYEn', 'o'); hold on;
xPlot = linspace(min(dataXEn(:)), max(dataXEn(:)), 100);
yPlot = hillEq(1, slop, 0, xPlot);
plot(10.^xPlot, yPlot, 'r'); xlabel('Relative Concentration'); ylabel('Norm.(\DeltaF/F)');
set(gca, 'XScale', 'log')
hold off;

saveas(gcf, fullfile('results', 'figures', 'Figure3A.fig'));

% plot Figure S4A
disp('----------REFINE PARAMETERS FOR SATURATED CURVES----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');
figure; 
for i = 1:length(ampVec)
    cMatrix(gfRC(i,1), gfRC(i,2)) = kdVec(i);
    aMatrix(gfRC(i,1), gfRC(i,2)) = ampVec(i);
    r2Matrix(gfRC(i,1), gfRC(i,2))= gfR2(i);
    
    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', odorList{gfRC(i,1)}, ...
        ORNList{gfRC(i, 2)}, ampVec(i), slop, kdVec(i), gfR2(i));
    
    xx = squeeze(conc(gfRC(i), ftCol(i), :));
    yy = squeeze(dff(ftRow(i), ftCol(i), :));
    sem = squeeze(dffSEM(ftRow(i), ftCol(i), :));


    subplot(6, 6, i);
    xP = linspace(min(log10(xx)), max(log10(xx)), 50);
    yP = hillEq( ampVec(i), slop, kdVec(i), xP);

    errorbar(log10(xx), yy, sem, 'ok', 'MarkerSize', 4); hold on;
    plot(xP, yP, 'k'); xlabel('log_{10}(c)'); ylabel('\DeltaF/F');
    title([odorList{gfRC(i, 1)}, ' / ', ORNList{gfRC(i, 2)}]);

end
set(gcf, 'Position', [63 1 1500 1000])

saveas(gcf, fullfile('results', 'figures', 'FigureS4A'));

%% Following code generate estimation of parameters for non-satruated curves
% The results will further be refined using MLE
% fit with fixed slop
ft2 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Robust = 'Bisquare';
opts2.Lower = [0 slop -11]; % setup the range and initial value of the variable
opts2.Upper = [max(dff(:))*1.5 slop 0];
opts2.StartPoint = [4 slop -7];

disp('----------FIT WITH KNOWING SLOP----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftCol, ftRow] = find(maskPlot' == 2);

for i = 1:length(ftRow)
    xx = squeeze(conc(ftRow(i), ftCol(i), :));
    yy = squeeze(dff(ftRow(i), ftCol(i), :));
    
    [fitresult, gof] = fit(log10(xx), yy, ft2, opts2);   %fit
    
    rSq = gof.rsquare; coeff = coeffvalues(fitresult);
    
    aMatrix(ftRow(i), ftCol(i)) = coeff(1);
    cMatrix(ftRow(i), ftCol(i)) = coeff(3);
    r2Matrix(ftRow(i), ftCol(i))= rSq;
       
    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', odorList{ftRow(i)}, ...
        ORNList{ftCol(i)}, coeff(1), coeff(2), coeff(3), rSq);
        
	if plotFlag == 1
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([odorList{ftRow(i)}, '/', ORNList{ftCol(i)}, '/', num2str(coeff(1)), '/', num2str(coeff(3))]);
    end
end

% Estimate the maximum values for each Odor
maxAofOdor = zeros(rowTotal, 1); % define parameter to store maximum amplitude 

disp('----------FIND AMPLITUDE OF EACH ODOR FOR WEAK RESPONSE CURVE FITTING----------');
fprintf('%30s\t%-5s\t\n', 'Odor Name', 'Amplitude');

 
for i = 1:rowTotal
    colIdx = ~isnan(aMatrix(i, :));
    maxSeq = aMatrix(i, colIdx); %find the nun-zero elements of each row (odor), from the saturated dataset

    maxAofOdor(i) = mean(maxSeq(~isoutlier(maxSeq)));

	fprintf('%30s\t%.2f\t\n', odorList{i}, maxAofOdor(i));
end

% fit with fixed slop and amplitude
ft3 = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Robust = 'Bisquare';

disp('----------FIT WITH KNOWING SLOP AND AMPLITUDE----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

[ftRow, ftCol] = find(maskPlot == 3);

for i = 1:length(ftRow)
    xx = squeeze(conc(ftRow(i), ftCol(i), :));
    yy = squeeze(dff(ftRow(i), ftCol(i), :));

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
       
    fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', odorList{ftRow(i)}, ...
        ORNList{ftCol(i)}, coeff(1), coeff(2), coeff(3), rSq);
        
	if plotFlag == 1
        xP = linspace(min(log10(xx)), max(log10(xx)), 50);
        yP = hillEq(coeff(1), coeff(2), coeff(3), xP);

        figure; 
        plot(log10(xx), yy, 'ok'); hold on;
        plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
        title([odorList{ftRow(i)}, ORNList{ftCol(i)}]);
    end
end

% show the failed fittings
gofThld = 0.4;

disp('----------FAILED FITTINGS----------');
fprintf('%25s\t%-10s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'Odor', 'ORN', 'Amp', 'Slop', 'EC50', 'R^2');

plotFlag = 0;
for i = 1: rowTotal
    for j = 1: colTotal
        if r2Matrix(i, j) < gofThld || cMatrix(i, j)  ==0
            
            maskPlot(i, j) = 4;

            fprintf('%25s\t%-5s\t%.2f\t%.2f\t%.2f\t%.2f\n', odorList{i}, ...
                ORNList{j}, aMatrix(i,j), slop, cMatrix(i,j), r2Matrix(i,j));

            if plotFlag == 1
                xx = squeeze(conc(i, j, :));
                yy = squeeze(dff(i, j, :));
    
                xP = linspace(min(log10(xx)), max(log10(xx)), 50);
                yP = hillEq( aMatrix(i,j), slop, cMatrix(i,j), xP);

                figure; 
                plot(log10(xx), yy, 'ok'); hold on;
                plot(xP, yP, 'r'); xlabel('log10(c)'); ylabel('\DeltaF/F');
                title([odorList{i}, ORNList{j}]);
            end
            
        end
    end
end

% save the current fitting results
save(fullfile('results','initialParameters.mat'), 'aMatrix', 'cMatrix', 'r2Matrix', 'slop', 'maskPlot');