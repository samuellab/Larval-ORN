%% input info
input.matFiles = {'13a_data.mat'; '22c_data.mat'; '42a_data.mat'}; %file names
input.varNames = {'sigMat'; 'pp'; 'pp'};  %variable names
input.ORNs = {'Or13a'; 'Or22c'; 'Or42a'}; %ORN name
input.odors = {'6-methyl-5-hepten-2-ol', '3-octanol'; 'methyl salicylate', 'anisole'; ...
    '4-hexen-3-one', '3-pentanol'}; %odor name;
input.concList = {[3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5; 10^-4], ...
    [10^-8; 3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5]; ...
    [10^-11; 3.16*10^-11; 10^-10; 3.16*10^-10; 10^-9; 3.16*10^-9; 10^-8; 3.16*10^-8], ...
    [3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5; 10^-4]; ...
    [10^-8; 3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5], ...
    [10^-8; 3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5]}; % concentation list

% other settings
cColor =[0 0.4470 0.7410; 0.85 0.325 0.0980];

%% output info
results.rMean = cell(3,2);
results.rVar = cell(3,2);
results.rSEM = cell(3,2);

results.fit = cell(3,2);

results.fitIndiv = cell(3, 2);
results.cmpPValue = cell(3, 1);

%%
for  ff = 1:2
    
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





