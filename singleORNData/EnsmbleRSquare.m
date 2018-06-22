function rSquareEn = EnsmbleRSquare(dataX, dataY, coefEn, funForm)
% get the goodness of fit, R-square, of the data ensemble

% calculate the predicted Y
[trialNum, inputLen] = size(dataX);
yPre = zeros(trialNum, inputLen);
for i = 1 : trialNum
    param = coefEn(i, :);
    xIn = dataX(i, :);
    yPre(i, :) = funForm(param(1), param(2), param(3), xIn);
end

% calculate sst, ssr and r-square
yVec    = dataY(:);
yPreVec = yPre(:);

yBar = mean(yVec);

sst = sum((yVec-yBar).^2);
ssr = sum((yVec-yPreVec).^2);

rSquareEn = 1-ssr/sst;

% plot the predict-vs-actual data
figure; 
plot(yVec, yPreVec, 'o'); hold on;

xPlot = linspace(min(yVec), max(yVec), 100);
plot(xPlot, xPlot, 'r'); hold off;
title(['{R^2} = ', num2str(rSquareEn)]);

end