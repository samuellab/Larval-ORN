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

%% other settings
% color map for curves
cColor = [252,141,89;...
    240,59,32; ...
    44,127,184; ...
    127,205,187]/255;

% setup fitting function and method
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -11]; % setup the range and initial value of the variable
opts.Upper = [10 15 -3];
opts.StartPoint = [5 5 -9];

%% averaged data
%calculate the mean and sem
results.rMean = cellfun(@mean, input.dff, 'UniformOutput', false);
results.rSEM  = cellfun(@(x) std(x)/size(x, 1), input.dff, 'UniformOutput', false);

% plot the curves in one figure
figure; hold on; legText = {};
for i = 1:length(input.ORNs)
    for j = 1:length(input.odors(1, :))
        errorbar(input.concList{i,j}, results.rMean{i, j}, results.rSEM{i, j}, ...
            '-o', 'Color', cColor(j+(i-1)*length(input.odors(1, :)),:)); 
        legText{end+1} = [input.ORNs{i}, ' ',input.odors{i, j}];
    end    
end
% format the figure 
legend(legText, 'Location',  'northwest'); axis tight; hold off;
xlabel('Concentration'); ylabel('\DeltaF/F'); set(gca,'XScale','log' );
yTickMax = round(max(reshape(cell2mat(results.rMean), [], 1)));
yticks(1:yTickMax); yticklabels(num2str((1:yTickMax)')); 
xticks(logspace(-11, -4, 8));
xticklabels({'10^{11}', '10^{10}', '10^{9}', '10^{8}', '10^{7}', ...
    '10^{6}', '10^{5}', '10^{4}',});
set(gcf, 'Position', [100, 100, 700, 420]); movegui(gcf, 'northwest');

% fit the averaged curuves and print
disp('----------FIT AVERAGED CURVE:----------');
fprintf('%5s\t%-20s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'ORN', 'Odor', 'Amp', 'Slop', 'EC_{50}', 'R^2');
for i = 1:length(input.ORNs)
    for j = 1:length(input.odors(1, :))
        [fitresult, gof] = fit(log10((input.concList{i,j})'), (results.rMean{i, j})', ft, opts);
        results.fitCoeff{i, j} = coeffvalues(fitresult);    
        results.fitR2{i, j} = gof.rsquare; 
        fprintf('%5s\t%-20s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            input.ORNs{i}, input.odors{i, j}, results.fitCoeff{i, j}(1),...
            results.fitCoeff{i, j}(2), results.fitCoeff{i, j}(3), results.fitR2{i, j});
    end
end

%% show each individual curve

% plot each individual curve on top of each other, seperate two ORNs
for i = 1:length(input.ORNs)
    figure; hold on;  title(input.ORNs{i});
    for j = 1:length(input.odors(1, :))
        plot(input.concList{i,j}, input.dff{i,j}, '-o', 'Color', ...
            cColor(j+(i-1)*length(input.odors(1, :)),:));
    end 
    
    axis tight; hold off;
    xlabel('Concentration'); ylabel('\DeltaF/F'); set(gca,'XScale','log' );    
end

%% fit each individual curve and display the parameters
disp('----------FIT INDIVIDUAL CURVE:----------');
fprintf('%5s\t%-20s\t%-5s\t%-5s\t%-5s\t%-5s\t%-5s\t\n', 'ORN', 'Odor', 'Trial', 'Amp', 'Slop', 'EC_{50}', 'R^2');
for i = 1:length(input.ORNs)
    for j = 1:length(input.odors(1, :))
        yMat = input.dff{i, j}; [trialNum, ~] = size(yMat);
        for k = 1:trialNum
            [fitresult, gof] = fit(log10((input.concList{i,j})'), (yMat(k,:))', ft, opts);
            results.fitCoeffIdv{i, j}(k,:) = coeffvalues(fitresult);    
            results.fitR2Idv{i, j}(k) = gof.rsquare; 
            fprintf('%5s\t%-20s\t%5s\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
                input.ORNs{i}, input.odors{i, j}, input.expID{i, j}{k}, ...
                results.fitCoeffIdv{i, j}(k,1), results.fitCoeffIdv{i, j}(k,2), ...
                results.fitCoeffIdv{i, j}(k,3), results.fitR2Idv{i, j}(k));
        end
    end
end

%% compare the parameters
f1 = figure; hold on; title('Amplitude');
f2 = figure; hold on; title('EC_{50}');
f3 = figure; hold on; title('Hill Coeff');
coeffPool = cell2mat(results.fitCoeffIdv);
ampPool = coeffPool(:, [1,4]); 
hcPool = coeffPool(:, [2,5]); 
kdPool = coeffPool(:, [3,6]); 
for i = 1:length(input.ORNs)
    for k = 1:length(input.expID{i, 1})
        myIndex = (i-1) * length(input.expID{1, 1}) + k;
        for j = 1:length(input.odors(1, :))
            figure(f1); plot((i-1)*length(input.odors(1, :)) + j, ampPool(myIndex, j), ...
                'o', 'Color', cColor(j+(i-1)*length(input.odors(1, :)),:));
            figure(f2); plot((i-1)*length(input.odors(1, :)) + j, hcPool(myIndex, j), ...
                'o', 'Color', cColor(j+(i-1)*length(input.odors(1, :)),:));
            figure(f3); plot((i-1)*length(input.odors(1, :)) + j, kdPool(myIndex, j), ...
                'o', 'Color', cColor(j+(i-1)*length(input.odors(1, :)),:));
        end
        figure(f1); plot([1 2]+(i-1)*2, ampPool(myIndex, :), 'k');
        figure(f2); plot([1 2]+(i-1)*2, hcPool(myIndex, :), 'k');
        figure(f3); plot([1 2]+(i-1)*2, kdPool(myIndex, :), 'k');
    end
end

figure(f1); xticks(1:i*j); hold off 
xticklabels({input.odors{1,1}, input.odors{1,2}, input.odors{2,1}, input.odors{2,2}});
xtickangle(-45);    ylabel('\DeltaF/F'); axis([0.5 4.5 floor(min(ampPool(:))) ceil(max(ampPool(:)))]);
set(gcf, 'Position', [100, 100, 350, 420]); movegui(gcf, 'north');

figure(f2); xticks(1:i*j); hold off 
xticklabels({input.odors{1,1}, input.odors{1,2}, input.odors{2,1}, input.odors{2,2}});
xtickangle(-45);    ylabel('\DeltaF/F'); axis([0.5 4.5 floor(min(hcPool(:))) ceil(max(hcPool(:)))]);
set(gcf, 'Position', [100, 100, 350, 420]); movegui(gcf, 'north');

figure(f3);xticks(1:i*j); hold off 
xticklabels({input.odors{1,1}, input.odors{1,2}, input.odors{2,1}, input.odors{2,2}});
xtickangle(-45);    ylabel('\DeltaF/F'); 
axis([0.5 4.5 floor(min(kdPool(:))) ceil(max(kdPool(:)))]);
set(gcf, 'Position', [100, 100, 350, 420]); movegui(gcf, 'north');

%% shift the EC_50, plot data
dffNorm = cellfun(@(x, y) x./y(:, 1), input.dff, results.fitCoeffIdv, 'UniformOutput', false);
concShift = cellfun(@(x, y) repmat(log10(x), length(y(:,3)), 1) - repmat(y(:, 3), 1, length(x)), ...
    input.concList, results.fitCoeffIdv, 'UniformOutput', false);

xData = []; yData = [];
for i = 1:length(input.ORNs)
    for j = 1:length(input.odors(1, :))
        xData = [xData; concShift{i, j}];
        yData = [yData; dffNorm{i, j}];
    end
end

figure; 
plot(xData', yData', 'ok'); hold off





























%%
for  ff = 2:2
    
    load(input.matFiles{ff}, input.varNames{ff});   % load files
    if  strcmp(input.varNames{ff}, 'pp')
        dataPool = pp';
    elseif strcmp(input.varNames{ff}, 'sigMat')
        dataPool = sigMat;
    else
        error('Not listed input variable.');
    end
    
    dataOdor1 = dataPool(:, 1:2:end-1); dataOdor2 = dataPool(:, 2:2:end);
    
    input.data{ff, 1} = dataOdor1;  input.data{ff, 2} = dataOdor2;
    
    [trialNum, concNum] = size(dataOdor1);

    

    %% fit on individual curve
    coefMat = zeros(trialNum, 6 );      

    px1 = linspace (log10(input.concList{ff,1}(1)), log10(input.concList{ff,1}(end)));
    px2 = linspace (log10(input.concList{ff,2}(1)), log10(input.concList{ff,2}(end)));

    figure; title(input.ORNs{ff});
    for i = 1:trialNum
        %fit
        [fitresult1, ~] = fit(log10(input.concList{ff,1}), dataOdor1(i, :)', ft, opts);
        [fitresult2, ~] = fit(log10(input.concList{ff,2}), dataOdor2(i, :)', ft, opts);

        %paramter
        coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
        coefMat(i, 1:3) = coef1; coefMat(i, 4:6) = coef2;

        py1 = hillEq(coef1(1), coef1(2), coef1(3), px1);
        py2 = hillEq(coef2(1), coef2(2), coef2(3), px2);

        subplot( ceil(trialNum/2), 2,  i);

        plot(input.concList{ff,1}, dataOdor1(i, :), 'o', 'Color', cColor(1,:));
        hold on;
        plot(10.^px1, py1,   'Color', cColor(1,:));

        plot(input.concList{ff,2}, dataOdor2(i, :), 'o', 'Color', cColor(2,:));
        plot(10.^px2, py2,   'Color', cColor(2,:));    

        xlabel('Concentration'); ylabel('\DeltaF/F'); 
        set(gca,'XScale','log' );

    %     legend({[odor1, ' ', num2str(coef1)], [odor2, ' ', num2str(coef2)]}, 'Location',  'northwest');
    %     title(['ORN trial=', num2str(i) ]);
        hold off;
    end
    
    % save the fitting results
    results.fitIndiv{ff, 1} = coefMat(:, 1:3);
    results.fitIndiv{ff, 2} = coefMat(:, 4:6);

    %% one-on-one compariation
    figure;
    cmpInd = 1:6;   cmpIndM= repmat(cmpInd, [trialNum, 1]);
    h = zeros(1, 3); p = h; % t-test score
    alpha = [1 1 -1];
    for i = 1:3
        [h(i), p(i)] = ttest(coefMat(:, i), coefMat(:, 3+i));
        plot(cmpIndM(:, 2*i-1:2*i)', alpha(i)*coefMat(:, [i, 3+i])', 'o-', 'color', 'k'); hold on;
    end
    xticks([1.5 3.5 5.5]);
    xticklabels({['Amp, h=', num2str(h(1))], ['Hill Coeff., h=', num2str(h(2))], ...
        ['EC_{50}, h=', num2str(h(3))]});
    title([input.ORNs{ff}, '-', input.odors{ff, 1}, ' vs ', input.odors{ff, 2}]);
    
    % save the comparision results
    results.cmpPValue{ff} =[h, p];

    %% get the averaged data

    data1Mean = mean(dataOdor1, 1);
    data1SEM = std(dataOdor1, 1)/sqrt(trialNum);
    data1Var = var(dataOdor1, 1);

    data2Mean = mean(dataOdor2, 1);
    data2SEM = std(dataOdor2, 1)/sqrt(trialNum);
    data2Var = var(dataOdor2, 1);
    
    results.rMean{ff, 1} = data1Mean;   results.rMean{ff, 2} = data2Mean;
    results.rSEM{ff, 1} = data1SEM;     results.rSEM{ff, 2} = data2SEM;
    results.rVar{ff, 1} = data1Var;     results.rVar{ff, 2} = data2Var;   
    
    %% fit the curve
    [fitresult1, gof1] = fit(log10(input.concList{ff,1}), data1Mean', ft, opts);
    [fitresult2, gof2] = fit(log10(input.concList{ff,2}), data2Mean', ft, opts);

    coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
    results.fit{ff, 1} = coef1;         results.fit{ff, 2} = coef2;
    
    py1 = hillEq(coef1(1), coef1(2), coef1(3), px1);
    py2 = hillEq(coef2(1), coef2(2), coef2(3), px2);

    %% plot the dose-response data
    figure; 
    errorbar(input.concList{ff,1}, data1Mean, data1SEM, 'o', 'Color', cColor(1,:)); hold on;
    plot(10.^px1, py1,   'Color', cColor(1,:));

    errorbar(input.concList{ff,2}, data2Mean, data2SEM, 'o', 'Color', cColor(2,:));
    plot(10.^px2, py2,   'Color', cColor(2,:));    

    legend({input.odors{ff,1}, input.odors{ff,2}}, 'Location',  'northwest');

    title([input.ORNs{ff} , 'N=', num2str(trialNum)])
    xlabel('Concentration'); ylabel('\DeltaF/F');
    set(gca,'XScale','log' );
    hold off;
end

% %% plot the var-vs-mean 
% or13aMean = [results.rMean{1,1}, results.rMean{1,2}];
% or13aVar  = [results.rVar{1,1}, results.rVar{1,2}];
% figure; 
% plot(or13aMean(:), or13aVar(:), 'ro');  hold on
% 
% or22cMean = [results.rMean{2,1}, results.rMean{2,2}];
% or22cVar  = [results.rVar{2,1}, results.rVar{2,2}];
% plot(or22cMean(:), or22cVar(:), 'bo'); 
% 
% or42aMean = [results.rMean{3,1}, results.rMean{3,2}];
% or42aVar  = [results.rVar{3,1}, results.rVar{3,2}];
% plot(or42aMean(:), or42aVar(:), 'ko'); hold off;

% %% plot the distribution of amp
% fitParmPool = cell2mat(results.fitIndiv);
% aPool = [fitParmPool(:, 1);  fitParmPool(:, 4)];
% figure; histogram(aPool); title('Amp');
% hPool = [fitParmPool(:, 2);  fitParmPool(:, 5)];
% figure; histogram(hPool); title('Hill Coeff');
% kdPool= [fitParmPool(:, 3);  fitParmPool(:, 6)];
% figure; histogram(kdPool); title('EC_{50}');

% %%
% fitParm = cell2mat(results.fitIndiv);
% ampPool = fitParm(:, [1, 4]);
% slopPool =fitParm(:, [2, 4]); 
% 
% kdPool1 = fitParm(1:12, 3);  
% kdPool2 = fitParm(1:12, 6); 
% 
% kdPool3 = fitParm(13:22, 3); 
% kdPool4 = fitParm(13:22, 6); 
% 
% kdPool5 = fitParm(23:33, 3); 
% kdPool6 = fitParm(23:33, 6); 
% 
% 
% ampMean = mean(ampPool(:))
% ampVar =  var(ampPool(:))
% 
% slopPool = slopPool(:);
% slopPool(slopPool>10) = [];
% slopMean = mean(slopPool(:))
% slopVar =  var(slopPool(:))
% 
% kdPoolCent = [kdPool1 - mean(kdPool1); kdPool2 - mean(kdPool2); kdPool3 - mean(kdPool3); ...
%     kdPool4 - mean(kdPool4); kdPool5 - mean(kdPool5); kdPool6 - mean(kdPool6)] ;
% kdVar = var(kdPoolCent)

%% discuss 13a's fitting results
% Q1, for the same neuron, is the measured \Delta{Kd} more reliable?
kd1 = results.fitIndiv{1,1}(:, 3);  kd2 = results.fitIndiv{1,2}(:, 3);
dkd = abs(kd1-kd2);

disp([input.ORNs{1}, ' responds to ', input.odors{1,1},' = ', ...
    num2str(mean(kd1)), ' +- ', num2str(std(kd1))]);
disp([input.ORNs{1}, ' responds to ', input.odors{1,2},' = ', ...
    num2str(mean(kd2)), ' +- ', num2str(std(kd2))]);
disp(['The difference of EC50 for each neuron: ', num2str(mean(dkd)), ...
    ' +- ', num2str(std(dkd))]);

% Q2, check if bilateral has less variation?
% calculate \Delta{Kd} between the two bilateral neurons
kd = [kd1, kd2];

kdL = kd(1:2:11, :); kdR = kd(2:2:12, :); 
dkdLR = kdL - kdR;
disp(['Bilateral difference of EC50: ', num2str(mean(dkdLR(:))), ' +- ', num2str(std(dkdLR(:)))]);

% Q3, is the amplitude more reliable for the same neuron? 
% compare to the fitting error to the difference. 

slop1 = results.fitIndiv{1,1}(:, 2);  slop2 = results.fitIndiv{1,2}(:, 2);
slopAll = [slop1, slop2];
mean(slopAll(:))
std(slopAll(:))





