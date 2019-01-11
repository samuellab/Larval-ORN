function [coeffs, gof] = Fit_DR(x, y, fitParam, plotFlag, figParam)
%Fit_DR(x, y, figName, plotFlag)
%  Fit the slop of dose-response curve using a logistic function.
%
%  Inputs:
%      x : independent variable data of the logistic function
%  	   y : dependent variable data of the logistic function
%      plotFlag: flag, 1 plot fitted curve; 0, do not plot
%      figHandle: the handle of figure to plot on
%      figName: str, figure name if need to plot plotFlag = 1
%  Outputs:
%      coeffs : vector, values of the fitted parameters
%      gof : structure, infor of goodness-of fit.

%% fit:  
[xData, yData] = prepareCurveData( x, y );

% Set up the equation and options.
ft = fittype( 'a/(1+ exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';

opts.Lower = [-1 -1 -10]; % setup the range and initial value of the variable
opts.Upper = [10 10 1];
opts.StartPoint = [4 0.1 -6];

if ~isempty(fitParam) % if there is fitting parameter passed in
    for i = 1:length(fitParam)
        p = fitParam(i);
        if ~isnan(p)  % if there is non NaN elements
            opts.Lower(i) = p;
            opts.Upper(i) = p;
            opts.StartPoint(i) = p;
        end
    end
end

% Fit the parameters in the equation
[fitresult, gof] = fit( xData, yData, ft, opts );
coeffs = coeffvalues(fitresult);

%% display 
if plotFlag==1  
	figure(figParam.handle);  %claim the figure to plot on

    hLine = plot(fitresult );   %plot the curve
    hLine.Color = figParam.lineColor;   %set the color
    set(hLine, 'DisplayName', figParam.lineName);   %set the name of the curve

    legend('off'); legend('show', 'Location', 'northeastoutside'); %control the display of legends
end

end