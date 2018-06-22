%% input info
input.matFiles = '22c_data.mat';
input.varNames = 'pp';
input.ORNs = 'Or22ca';
input.odors = {'methyl salicylate', 'anisole'};
    
input.concList = {[10^-11; 3.16*10^-11; 10^-10; 3.16*10^-10; 10^-9; 3.16*10^-9; 10^-8; 3.16*10^-8], ...
    [3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5; 10^-4];}; % concentation list

% other settings
cColor =[0 0.4470 0.7410; 0.85 0.325 0.0980];

%% output info
results.fitCoef = cell(1, 2);
results.fitConfintWidth = cell(1, 2);

results.fitCoefClean = cell(1, 2);
results.fitConfintWidthClean = cell(1, 2);

results.fitCoefGroup = cell(2, 2);
results.fitConfintWidthGroup = cell(2, 2);

%% load data and setup fitting
load(input.matFiles, input.varNames);   % load files

if  strcmp(input.varNames, 'pp')
    dataPool = pp';
elseif strcmp(input.varNames, 'sigMat')
    dataPool = sigMat;
else
    error('Not listed input variable.');
end

dataOdor1 = dataPool(:, 1:2:end-1); dataOdor2 = dataPool(:, 2:2:end);

input.data{1, 1} = dataOdor1;  input.data{1, 2} = dataOdor2;

[trialNum, concNum] = size(dataOdor1);

resp = cell2mat(input.data);
dataY = [resp(:, 1:concNum); resp(:, concNum+1:2*concNum)];

conc = cell2mat(input.concList);
conc = conc';
dataX = [repmat(conc(1,:), [trialNum,1]); repmat(conc(2,:), [trialNum,1])];
dataX = log10(dataX);

% setup fitting function and method
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -12]; % setup the range and initial value of the variable
opts.Upper = [10 15 1];
opts.StartPoint = [4 5 -6.5];

%% fit on individual curve

px1 = linspace (log10(input.concList{1,1}(1)), log10(input.concList{1,1}(end)));
px2 = linspace (log10(input.concList{1,2}(1)), log10(input.concList{1,2}(end)));

figure; title(input.ORNs);
sse = 0;
for i = 1:trialNum
    %fit
    [fitresult1, gof1] = fit(log10(input.concList{1,1}), dataOdor1(i, :)', ft, opts);
    [fitresult2, gof2] = fit(log10(input.concList{1,2}), dataOdor2(i, :)', ft, opts);

    %paramter
    coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
    results.fitCoef{1}(i, :) = coef1; results.fitCoef{2}(i, :) = coef2;
    
    ci1 = confint(fitresult1); ciw1 = (ci1(2, :) - ci1(1, :))/2; 
    ci2 = confint(fitresult2); ciw2 = (ci2(2, :) - ci2(1, :))/2; 
    results.fitConfintWidth{1}(i, :) = ciw1; results.fitConfintWidth{2}(i, :) = ciw2;
    
    sse = sse + gof1.sse + gof2.sse;

    py1 = hillEq(coef1(1), coef1(2), coef1(3), px1);
    py2 = hillEq(coef2(1), coef2(2), coef2(3), px2);

    subplot( ceil(trialNum/2), 2,  i);

    plot(input.concList{1,1}, dataOdor1(i, :), 'o', 'Color', cColor(1,:));
    hold on;
    plot(10.^px1, py1,   'Color', cColor(1,:));

    plot(input.concList{1,2}, dataOdor2(i, :), 'o', 'Color', cColor(2,:));
    plot(10.^px2, py2,   'Color', cColor(2,:));    

    xlabel('Concentration'); ylabel('\DeltaF/F'); 
    set(gca,'XScale','log' );

    hold off;
end
rSqIndFit = 1- sse/sum((dataY(:) - mean(dataY(:))).^2);
disp(['R-Square of fitting each individual curve = ',num2str(rSqIndFit)]);

%% display the fitting parameters
lrStr = {'R', 'L'};
disp('----------FITTING RESULTS w/s ALL DATA:----------');
fprintf('%10s\t%-5s\t%-5s\t%-10s\t%-10s\t%-10s\t\n', 'Trial#', 'L/R', 'Odor', 'Amp', 'Slop', 'EC_{50}');

for i = 1 : trialNum
    trilNum = ceil(i/2);
    lrMark = rem(i, 2) + 1;
    
    for j = 1:2
       a = results.fitCoef{j}(i, 1);
       da = results.fitConfintWidth{j}(i, 1);
       b = results.fitCoef{j}(i, 2);
       db = results.fitConfintWidth{j}(i, 2);
       c = results.fitCoef{j}(i, 3);
       dc = results.fitConfintWidth{j}(i, 3);
 
       fprintf('%10d\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n',trilNum, lrStr{lrMark}, j, a, da, b, db, c, dc);
       
    end
end

cMatTemp = cell2mat(results.fitCoef);
disp('------------------------------------------------------');
fprintf('%10s\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n', 'Average', '-', 1,  ...
    mean(cMatTemp(:, 1)), std(cMatTemp(:, 1))/sqrt(trialNum), ...
    mean(cMatTemp(:, 2)), std(cMatTemp(:, 2))/sqrt(trialNum), ...
    mean(cMatTemp(:, 3)), std(cMatTemp(:, 3))/sqrt(trialNum));
fprintf('%10s\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n', 'Average', '-', 2,  ...
    mean(cMatTemp(:, 4)), std(cMatTemp(:, 4))/sqrt(trialNum), ...
    mean(cMatTemp(:, 5)), std(cMatTemp(:, 5))/sqrt(trialNum), ...
    mean(cMatTemp(:, 6)), std(cMatTemp(:, 6))/sqrt(trialNum));

%%  remove 22c's decayed data
input.satMaskO1 = [...
    1 1 1 1 1 1 1 1; ...    %1L
    1 1 1 1 1 1 1 1; ...    %1R
    1 1 1 1 1 1 1 1; ...    %2L
    1 1 1 1 1 1 1 1; ...    %2R
    1 1 1 1 1 1 1 1; ...    %3L
    1 1 1 1 1 1 1 1; ...    %3R
    1 1 1 1 1 1 1 1; ...    %4L
    1 1 1 1 1 1 1 1; ...    %4R
    1 1 1 1 1 1 1 1; ...    %5L
    1 1 1 1 1 1 1 1];       %5R

input.satMaskO2 = [...
    1 1 1 1 1 1 1 0; ...    %1L
    1 1 1 1 1 1 1 0; ...    %1R
    1 1 1 1 1 1 1 0; ...    %2L
    1 1 1 1 1 1 1 0; ...    %2R
    1 1 1 1 1 1 1 0; ...    %3L
    1 1 1 1 1 1 1 0; ...    %3R
    1 1 1 1 1 1 1 0; ...    %4L
    1 1 1 1 1 1 0 0; ...    %4R
    1 1 1 1 1 1 1 0; ...    %5L
    1 1 1 1 1 1 0 0];       %5R

% fit on individual curve with mask
figure; title(['Cleaned ',input.ORNs]);

for i = 1:trialNum
    % get the index defined by the mask
    usedList1 = find(input.satMaskO1(i, :));
    usedList2 = find(input.satMaskO2(i, :));
    
    in1 = log10(input.concList{1,1}(usedList1));
    in2 = log10(input.concList{1,2}(usedList2));
    
    out1 = dataOdor1(i, usedList1)'; 
    out2 = dataOdor2(i, usedList2)'; 
    
    [fitresult1, ~] = fit(in1, out1, ft, opts);
    [fitresult2, ~] = fit(in2, out2, ft, opts);

    %paramter
    coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
    results.fitCoefClean{1}(i, :) = coef1;
    results.fitCoefClean{2}(i, :) = coef2;
    
    ci1 = confint(fitresult1); ciw1 = (ci1(2, :) - ci1(1, :))/2; 
    ci2 = confint(fitresult2); ciw2 = (ci2(2, :) - ci2(1, :))/2; 
    results.fitConfintWidthClean{1}(i, :) = ciw1; results.fitConfintWidthClean{2}(i, :) = ciw2;

    py1 = hillEq(coef1(1), coef1(2), coef1(3), px1);
    py2 = hillEq(coef2(1), coef2(2), coef2(3), px2);

    subplot( ceil(trialNum/2), 2,  i);

    plot(10.^in1, out1, 'o', 'Color', cColor(1,:));
    hold on;
    plot(10.^px1, py1,   'Color', cColor(1,:));

    plot(10.^in2, out2, 'o', 'Color', cColor(2,:));
    plot(10.^px2, py2,   'Color', cColor(2,:));    

    xlabel('Concentration'); ylabel('\DeltaF/F'); 
    set(gca,'XScale','log' );

    hold off;
end

%%
disp('----------FITTING RESULTS w/s MASKED DATA:----------');
fprintf('%10s\t%-5s\t%-5s\t%-10s\t%-10s\t%-10s\t\n', 'Trial#', 'L/R', 'Odor', 'Amp', 'Slop', 'EC_{50}');

for i = 1 : trialNum
    trilNum = ceil(i/2);
    lrMark = rem(i, 2) + 1;
    
    for j = 1:2
       a = results.fitCoefClean{j}(i, 1);
       da = results.fitConfintWidthClean{j}(i, 1);
       b = results.fitCoefClean{j}(i, 2);
       db = results.fitConfintWidthClean{j}(i, 2);
       c = results.fitCoefClean{j}(i, 3);
       dc = results.fitConfintWidthClean{j}(i, 3);
 
       fprintf('%10d\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n',trilNum, lrStr{lrMark}, j, a, da, b, db, c, dc);
       
    end
end

cMatMaskTemp = cell2mat(results.fitCoefClean);
disp('------------------------------------------------------');
fprintf('%10s\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n', 'Average', '-', 1,  ...
    mean(cMatMaskTemp(:, 1)), std(cMatMaskTemp(:, 1))/sqrt(trialNum), ...
    mean(cMatMaskTemp(2:end, 2)), std(cMatMaskTemp(2:end, 2))/sqrt(trialNum-1), ...
    mean(cMatMaskTemp(:, 3)), std(cMatMaskTemp(:, 3))/sqrt(trialNum));
fprintf('%10s\t%-5s\t%-5d\t%.2f+-%.2f\t%.2f+-%.2f\t%.2f+-%.2f\t\n', 'Average', '-', 2,  ...
    mean(cMatMaskTemp(:, 4)), std(cMatMaskTemp(:, 4))/sqrt(trialNum), ...
    mean(cMatMaskTemp(:, 5)), std(cMatMaskTemp(:, 5))/sqrt(trialNum), ...
    mean(cMatMaskTemp(:, 6)), std(cMatMaskTemp(:, 6))/sqrt(trialNum));

%% side-by-side compare the fitted variables before and after cleanning
cMat = [cMatTemp(:, 1:3); cMatTemp(:, 4:6)];
cMatMask = [cMatMaskTemp(:, 1:3); cMatMaskTemp(:, 4:6)];

figure;	alpha = [1 1 -1];
for i = 1:3
    cmpInd = 1:6;   cmpIndM= repmat(cmpInd, [trialNum*2, 1]);

	plot(cmpIndM(:, 2*i-1:2*i)', alpha(i)*[cMat(:, i),  cMatMask(:, i)]', 'o-', 'color', 'k'); 
    hold on;
end

xticks([1.5 3.5 5.5]);  xticklabels({'Amp',  'Hill Coeff.', 'EC_{50}'});
title('All data vs Masked data');


%% compare the \Delta{mean} with averaged mean
ratio = cellfun(@(x, y) abs(x-y)./((x+y)/2), results.fitCoef, results.fitCoefClean, 'UniformOutput', false);

disp('----------PERCENTAGE CHANGE BEFORE & AFTER CLEANING----------');
fprintf('%10s\t%-5s\t%-5s\t%-10s\t%-10s\t%-10s\t\n', 'Trial#', 'L/R', 'Odor', 'Amp', 'Slop', 'EC_{50}');

for i = 1 : trialNum
    trilNum = ceil(i/2);        lrMark = rem(i, 2) + 1;
    for j = 1:2
       ra = ratio{j}(i, 1);     rb = ratio{j}(i, 2);    rc = ratio{j}(i, 3);   
       fprintf('%10d\t%-5s\t%-5d\t%.3f\t%.3f\t%.3f\t\n',trilNum, lrStr{lrMark}, j, ra, rb, rc);
    end
end

rMat = cell2mat(ratio);     %ratio matrix
rMat = [rMat(2:end, 1:3); rMat(:, 4:6)];

rMean = mean(rMat, 1);
disp('------------------------------------------------------');
fprintf('%10s\t%-5s\t%-5s\t%.3f\t%.3f\t%.3f\t\n', 'Average', '-', '-',  rMean(1), rMean(2), rMean(3));

%% undersample the data and check how the fittings change
g1 = 1:2:length(input.concList{1})-1;
g2 = 2:2:length(input.concList{1});
g3 = [1 2 4 6];
group = [g1; g2; g3];


for g = 1:3
    figure; title([input.ORNs, ' Group ', num2str(g)]);
    
    for i = 1:trialNum
        %fit
        [fitresult1, ~] = fit(log10(input.concList{1,1}(group(g, :))), dataOdor1(i, group(g, :))', ft, opts);
        [fitresult2, ~] = fit(log10(input.concList{1,2}(group(g, :))), dataOdor2(i, group(g, :))', ft, opts);

        %parameter
        coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
        results.fitCoefGroup{g, 1}(i, :) = coef1; results.fitCoefGroup{g, 2}(i, :) = coef2;

        ci1 = confint(fitresult1); ciw1 = (ci1(2, :) - ci1(1, :))/2; 
        ci2 = confint(fitresult2); ciw2 = (ci2(2, :) - ci2(1, :))/2; 
        results.fitConfintWidth{g, 1}(i, :) = ciw1; results.fitConfintWidth{g,2}(i, :) = ciw2;

        py1 = hillEq(coef1(1), coef1(2), coef1(3), px1);
        py2 = hillEq(coef2(1), coef2(2), coef2(3), px2);

        subplot( ceil(trialNum/2), 2,  i);

        plot(input.concList{1,1}(group(g, :)), dataOdor1(i, group(g, :)), 'o', 'Color', cColor(1,:));
        hold on;
        plot(10.^px1, py1,   'Color', cColor(1,:));

        plot(input.concList{1,2}(group(g, :)), dataOdor2(i, group(g, :)), 'o', 'Color', cColor(2,:));
        plot(10.^px2, py2,   'Color', cColor(2,:));    

        xlabel('Concentration'); ylabel('\DeltaF/F'); 
        set(gca,'XScale','log' );

        hold off;
    end
end

%%
disp('----------FITTING RESULTS UNDERSAMPLE:----------');
fprintf('%10s\t%-5s\t%-5s\t%-20s\t%-20s\t%-20s\t\n', 'Trial#', 'L/R', 'Odor', 'Amp', 'Slop', 'EC_{50}');

for i = 1 : trialNum
    trilNum = ceil(i/2);
    lrMark = rem(i, 2) + 1;
    
    for j = 1:2
       a = results.fitCoefClean{j}(i, 1);
       ag1 = results.fitCoefGroup{1, j}(i, 1);
       ag2 = results.fitCoefGroup{2, j}(i, 1);
       ag3 = results.fitCoefGroup{3, j}(i, 1);
       b = results.fitCoefClean{j}(i, 2);
       bg1 = results.fitCoefGroup{1, j}(i, 2);
       bg2 = results.fitCoefGroup{2, j}(i, 2);
       bg3 = results.fitCoefGroup{3, j}(i, 2);
       c = results.fitCoefClean{j}(i, 3);
       cg1 = results.fitCoefGroup{1, j}(i, 3);
       cg2 = results.fitCoefGroup{2, j}(i, 3);
       cg3 = results.fitCoefGroup{3, j}(i, 3);
 
       fprintf('%10d\t%-5s\t%-5d\t%.2f/%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f/%.2f\t%.2f/%.2f/%.2f/%.2f\t\n',trilNum, lrStr{lrMark}, j, a, ag1, ag2, ag3, b, bg1, bg2, bg3, c, cg1, cg2, cg3);
       
    end
end

%%
cTemp = cell2mat(results.fitCoefClean);
cgTemp = cell2mat(results.fitCoefGroup);

c = cTemp(: , [3, 6]);	cc = c(:);
cg= cgTemp(:, [3, 6]);

for i = 1:3
    cgg = cg(trialNum*(i-1)+1 : trialNum*i, :);
    
    ccgg = cgg(:);
    ddc(:, i)  = abs(cc - ccgg);
end
ddc(2:end,:)
mean(ddc(2:end,:), 1)
std(ddc(2:end,:),1)

% so, estimation of ec-50 is robust of the sampling

%% ensemble fitting of all data


[slop, ampVec, kdVec, dSlop, rSqaureSlop, rSquareEn] = EnsembleFit(dataX, dataY);

disp('Ensemble fitting of both odors');
disp(['slop = ', num2str(slop), ' +- ', num2str((dSlop(2)-dSlop(1))/2)]);
disp(['Kd_1 = ', num2str(mean(kdVec(1:trialNum))), ' +- ',  num2str(std(kdVec(1:trialNum))/sqrt(trialNum))]);
disp(['Kd_2 = ',num2str(mean(kdVec(1+trialNum:2*trialNum))), ' +- ',  num2str(std(kdVec(trialNum+1:2*trialNum))/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVec(1:trialNum))), ' +- ',  num2str(std(ampVec(1:trialNum))/sqrt(trialNum))]);
disp(['Amp_2 = ',num2str(mean(ampVec(1+trialNum:2*trialNum))), ' +- ',  num2str(std(ampVec(trialNum+1:2*trialNum))/sqrt(trialNum))]);
disp(['R-Square of the shared slop fit = ',num2str(rSqaureSlop)]);
disp(['R-Square of the ensemble fit = ',num2str(rSquareEn)]);

%% ensemble fitting on seperated odor data
% resp = cell2mat(input.data);
dataY1 = resp(:, 1:concNum);
dataY2 = resp(:, concNum+1: 2*concNum);

% conc = cell2mat(input.concList);
% conc = conc';
dataX1 = log10(repmat(conc(1,:), [trialNum,1])); 
dataX2 = log10(repmat(conc(2,:), [trialNum,1]));
 
[slop1, ampVec1, kdVec1, dSlop1, rSqaure2nd1, rSquareEn1] = EnsembleFit(dataX1, dataY1);
[slop2, ampVec2, kdVec2, dSlop2, rSqaure2nd2, rSquareEn2] = EnsembleFit(dataX2, dataY2);

disp('Ensemble fitting of odor 1');
disp(['slop_1 = ', num2str(slop1), ' +- ', num2str((dSlop1(2)-dSlop1(1))/2)]);
disp(['Kd_1 = ', num2str(mean(kdVec1)), ' +- ',  num2str(std(kdVec1)/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVec1)), ' +- ',  num2str(std(ampVec1)/sqrt(trialNum))]);
disp(['R-Square of the shared slop fit = ',num2str(rSqaure2nd1)]);
disp(['R-Square of the ensemble fit = ',num2str(rSquareEn1)]);

disp('Ensemble fitting of odor 2');
disp(['slop_2 = ', num2str(slop2), ' +- ', num2str((dSlop2(2)-dSlop2(1))/2)]);
disp(['Kd_1 = ', num2str(mean(kdVec2)), ' +- ',  num2str(std(kdVec2)/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVec2)), ' +- ',  num2str(std(ampVec2)/sqrt(trialNum))]);
disp(['R-Square of the shared slop fit = ',num2str(rSqaure2nd2)]);
disp(['R-Square of the ensemble fit = ',num2str(rSquareEn2)]);

%% split data into two parts, to mimic undersampled concentration
dataXG1 = dataX(:, 1:2:end-1);  dataXG2 = dataX(:, 2:2:end);
dataYG1 = dataY(:, 1:2:end-1);  dataYG2 = dataY(:, 2:2:end);

[slopG1, ampVecG1, kdVecG1, dSlopG1, rSqaure2ndG1, rSquareEnG1] = EnsembleFit(dataXG1, dataYG1);
[slopG2, ampVecG2, kdVecG2, dSlopG2, rSqaure2ndG2, rSquareEnG2] = EnsembleFit(dataXG2, dataYG2);

disp('Ensemble fitting of data group 1');
disp(['slop = ', num2str(slopG1), ' +- ', num2str((dSlopG1(2)-dSlopG1(1))/2)]);
disp(['Kd_1 = ', num2str(mean(kdVecG1(1:trialNum))), ' +- ',  num2str(std(kdVecG1(1:trialNum))/sqrt(trialNum))]);
disp(['Kd_2 = ',num2str(mean(kdVecG1(1+trialNum:2*trialNum))), ' +- ',  num2str(std(kdVecG1(1+trialNum:2*trialNum))/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVecG1(1:trialNum))), ' +- ',  num2str(std(ampVecG1(1:trialNum))/sqrt(trialNum))]);
disp(['Amp_2 = ',num2str(mean(ampVecG1(1+trialNum:2*trialNum))), ' +- ',  num2str(std(ampVecG1(1+trialNum:2*trialNum))/sqrt(trialNum))]);
disp(['R-Square of the shared slop fit = ',num2str(rSqaure2ndG1)]);
disp(['R-Square of the ensemble fit = ',num2str(rSquareEnG1)]);

disp('Ensemble fitting of data group 2');
disp(['slop = ', num2str(slopG2), ' +- ', num2str((dSlopG2(2)-dSlopG2(1))/2)]);
disp(['Kd_1 = ', num2str(mean(kdVecG2(1:trialNum))), ' +- ',  num2str(std(kdVecG2(1:trialNum))/sqrt(trialNum))]);
disp(['Kd_2 = ',num2str(mean(kdVecG2(1+trialNum:2*trialNum))), ' +- ',  num2str(std(kdVecG2(1+trialNum:2*trialNum))/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVecG2(1:trialNum))), ' +- ',  num2str(std(ampVecG2(1:12))/sqrt(trialNum))]);
disp(['Amp_2 = ',num2str(mean(ampVecG2(1+trialNum:2*trialNum))), ' +- ',  num2str(std(ampVecG2(1+trialNum:2*trialNum))/sqrt(trialNum))]);
disp(['R-Square of the shared slop fit = ',num2str(rSqaure2ndG2)]);
disp(['R-Square of the ensemble fit = ',num2str(rSquareEnG2)]);


%% use fminsearch to fit the parameters
% fit all data
[slopfs, ampVecfs, kdVecfs, rSquarefs] = EnsembleMiniSearch ...
    (dataX, dataY, hillEq, slop, ampVec, kdVec);

disp('--------------------fminsearch fitting:--------------------');
disp(['slop = ', num2str(slopfs)]);
disp(['Kd_1 = ', num2str(mean(kdVecfs(1:trialNum))), ' +- ',  num2str(std(kdVecfs(1:trialNum))/sqrt(trialNum))]);
disp(['Kd_2 = ',num2str(mean(kdVecfs(trialNum+1:2*trialNum))), ' +- ',  num2str(std(kdVecfs(trialNum+1:2*trialNum))/sqrt(trialNum))]);
disp(['Amp_1 = ', num2str(mean(ampVecfs(1:trialNum))), ' +- ',  num2str(std(ampVecfs(1:trialNum))/sqrt(trialNum))]);
disp(['Amp_2 = ',num2str(mean(ampVecfs(1+trialNum:2*trialNum))), ' +- ',  num2str(std(ampVecfs(1+trialNum:2*trialNum))/sqrt(trialNum))]);
disp(['R-Square = ',num2str(rSquarefs)]);

%% fit the undersampled data
[slopfsG1, ampVecfsG1, kdVecfsG1, rSquarefsG1] = EnsembleMiniSearch ...
    (dataXG1, dataYG1, hillEq, slopG1, ampVecG1, kdVecG1);

disp('--------------------fminsearch fitting group 1:--------------------');
disp(['slop = ', num2str(slopfsG1)]);
disp(['R-Square = ',num2str(rSquarefsG1)]);


[slopfsG2, ampVecfsG2, kdVecfsG2, rSquarefsG2] = EnsembleMiniSearch ...
    (dataXG2, dataYG2, hillEq, slopG2, ampVecG2, kdVecG2);

disp('--------------------fminsearch fitting group 2:--------------------');
disp(['slop = ', num2str(slopfsG2)]);
disp(['R-Square = ',num2str(rSquarefsG2)]);
