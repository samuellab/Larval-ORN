load('compiled_13a_data.mat', 'sigMat');

ORN = 'Or13a';

odor1 = '6-methyl-5-hepten-2-ol';
odor2 = '3-octanol';

concListOdor1 = [3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5; 10^-4];
concListOdor2 = [10^-8; 3.16*10^-8; 10^-7; 3.16*10^-7; 10^-6; 3.16*10^-6; 10^-5; 3.16*10^-5];

dataOdor1 = sigMat(:, 1:2:end-1);
dataOdor2 = sigMat(:, 2:2:end);

[trialNum, concNum] = size(dataOdor1);

%%
cColor =[0 0.4470 0.7410; 0.85 0.325 0.0980];
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -10]; % setup the range and initial value of the variable
opts.Upper = [10 10 1];
opts.StartPoint = [5 1 -6];
 
%% fit on individual curve
coefMat = zeros(trialNum, 6 );      

for i = 1:trialNum
    %fit
    [fitresult1, ~] = fit(log10(concListOdor1), dataOdor1(i, :)', ft, opts);
    [fitresult2, ~] = fit(log10(concListOdor2), dataOdor2(i, :)', ft, opts);
    
    %paramter
    coef1 = coeffvalues(fitresult1);    coef2 = coeffvalues(fitresult2);
    coefMat(i, 1:3) = coef1; coefMat(i, 4:6) = coef2;

%     figure; 
%     plot(concListOdor1, dataOdor1(i, :), 'o', 'Color', cColor(1,:));
%     hold on;
%     plot(concListOdor2, dataOdor2(i, :), 'o', 'Color', cColor(2,:));
% 
%     xlabel('Concentration'); ylabel('\DeltaF/F');
%     set(gca,'XScale','log' );
% 
%     legend({[odor1, ' ', num2str(coef1)], [odor2, ' ', num2str(coef2)]}, 'Location',  'northwest');
%     
%     title(['ORN trial=', num2str(i) ]);
%     hold off;
end




%% fit on the averaged data

data1Mean = mean(dataOdor1, 1);
data1SEM = std(dataOdor1, 1)/sqrt(trialNum);

data2Mean = mean(dataOdor2, 1);
data2SEM = std(dataOdor2, 1)/sqrt(trialNum);

%% plot the dose-response data
figure; 
errorbar(concListOdor1, data1Mean, data1SEM, 'Color', cColor(1,:));
hold on;
errorbar(concListOdor2, data2Mean, data2SEM, 'Color', cColor(2,:));

legend({odor1, odor2}, 'Location',  'northwest');

title(['ORN N=', num2str(trialNum)])
xlabel('Concentration'); ylabel('\DeltaF/F');
set(gca,'XScale','log' );
hold off;

%% fit the curve


[fitresult1, gof1] = fit(log10(concListOdor1), data1Mean', ft, opts);
[fitresult2, gof2] = fit(log10(concListOdor2), data2Mean', ft, opts);

% figure; plot(fitresult1);
% figure; plot(fitresult2);

%% compare the ymax
coef1 = coeffvalues(fitresult1);
coef2 = coeffvalues(fitresult2);

yMax1 = coef1(1); yMax2 = coef2(1); 

disp(['y_max1 = ', num2str(yMax1)]);
disp(['y_max2 = ', num2str(yMax2)]);


%%
