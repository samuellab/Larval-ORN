function [slop, ampVec, kdVec] = EnsembleFit(dataX, dataY)
  
[trialNum, xPoints] = size(dataX);

%% setup fitting function and method
hillEq = @(a, b, c, x)  a./(1+ exp(-b*(x-c)));
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.Lower = [-1 -1 -12]; % setup the range and initial value of the variable
opts.Upper = [10 30 1];
opts.StartPoint = [4 5 -6.5];

%% 1st round fit, get the Amp and EC_50 for each
coefR1 = zeros(trialNum, 3);
for i = 1:trialNum
    %fit
    [fitresult, ~] = fit(dataX(i, :)', dataY(i, :)', ft, opts);

    %paramter
    coefR1(i, :) = coeffvalues(fitresult);    
end

ampVec = coefR1(:, 1);

%% normalize amp and shift EC_50
dataYNorm = dataY./repmat(coefR1(:,1), [1, xPoints]);
dataXShift = dataX - repmat(coefR1(:,3), [1, xPoints]);

dataYNormVec = dataYNorm(:);
dataXShiftVec = dataXShift(:);

figure;
plot(dataXShiftVec, dataYNormVec, 'o'); hold on;

[fitresult, ~] = fit(dataXShiftVec , dataYNormVec , ft, opts);
coefR2 = coeffvalues(fitresult);
slop = coefR2(2);
xPlot = linspace(-2, 3, 100 );
yPlot =  hillEq(coefR2(1), coefR2(2), coefR2(3), xPlot);
plot(xPlot, yPlot);

%%
opts.Lower = [1 slop -12]; % setup the range and initial value of the variable
opts.Upper = [1 slop 1];
opts.StartPoint = [1 slop -6.5];

coefR3 = zeros(trialNum, 3);
for i = 1:trialNum
    %fit
    [fitresult, ~] = fit(dataX(i, :)', dataYNorm(i, :)', ft, opts);

    %paramter
    coefR3(i, :) = coeffvalues(fitresult);    
end

kdVec = coefR3(:, 3);
end
